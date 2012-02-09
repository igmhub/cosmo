// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/cosmo.h"
#include "likely/likely.h"
// the following are not part of the public API, so not included by likely.h
#include "likely/MinuitEngine.h"
#include "likely/EngineRegistry.h"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"

#include "boost/program_options.hpp"
#include "boost/bind.hpp"
#include "boost/ref.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "boost/format.hpp"
#include "boost/foreach.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

namespace lk = likely;
namespace po = boost::program_options;

class BaoFitPower {
public:
    BaoFitPower(cosmo::PowerSpectrumPtr fiducial, cosmo::PowerSpectrumPtr nowiggles)
    : _fiducial(fiducial), _nowiggles(nowiggles) {
        assert(fiducial);
        assert(nowiggles);
        setAmplitude(1);
        setScale(1);
        setSigma(0);
    }
    // Setter methods
    void setAmplitude(double value) { _amplitude = value; }
    void setScale(double value) { _scale = value; double tmp(value*value); _scale4 = tmp*tmp; }
    void setSigma(double value) { _sigma = value; _sigma2 = value*value; }
    // Returns the hybrid power k^3/(2pi^2) P(k) at the specified wavenumber k in Mpc/h.
    double operator()(double k) const {
        double ak(k/_scale), smooth(std::exp(-ak*ak*_sigma2/2));
        double fiducialPower = (*_fiducial)(ak), nowigglesPower = (*_nowiggles)(ak);
        return _scale4*(_amplitude*smooth*(fiducialPower - nowigglesPower) + nowigglesPower);
    }
private:
    double _amplitude, _scale, _scale4, _sigma, _sigma2;
    cosmo::PowerSpectrumPtr _fiducial, _nowiggles;
}; // BaoFitPower

typedef boost::shared_ptr<BaoFitPower> BaoFitPowerPtr;

class Binning {
public:
    Binning(int nBins, double lowEdge, double binSize)
    : _nBins(nBins), _lowEdge(lowEdge), _binSize(binSize) {
        assert(nBins > 0);
        assert(binSize > 0);
    }
    // Returns the bin index [0,nBins-1] or else -1.
    int getBinIndex(double value) const {
        int bin = std::floor((value - _lowEdge)/_binSize);
        assert(bin >= 0 && bin < _nBins);
        return bin;
    }
    // Returns the midpoint value of the specified bin.
    double getBinCenter(int index) const {
        assert(index >= 0 && index < _nBins);
        return _lowEdge + (index+0.5)*_binSize;
    }
    int getNBins() const { return _nBins; }
    double getLowEdge() const { return _lowEdge; }
    double getBinSize() const { return _binSize; }
private:
    int _nBins;
    double _lowEdge, _binSize;
}; // Binning

typedef boost::shared_ptr<const Binning> BinningPtr;

BinningPtr oversampleBinning(Binning const &other, int factor) {
    assert(factor > 0);
    BinningPtr bptr(new Binning(other.getNBins()*factor,other.getLowEdge(),other.getBinSize()/factor));
    return bptr;
}

class LyaData {
public:
    LyaData(BinningPtr logLambdaBinning, BinningPtr separationBinning, BinningPtr redshiftBinning,
    cosmo::AbsHomogeneousUniversePtr cosmology) : _cosmology(cosmology), _logLambdaBinning(logLambdaBinning),
    _separationBinning(separationBinning), _redshiftBinning(redshiftBinning)
    {
        assert(logLambdaBinning);
        assert(separationBinning);
        assert(redshiftBinning);
        assert(cosmology);
        _nsep = separationBinning->getNBins();
        _nz = redshiftBinning->getNBins();
        int nBinsTotal = logLambdaBinning->getNBins()*_nsep*_nz;
        _data.resize(nBinsTotal,0);
        _cov.resize(nBinsTotal,0);
        _r3d.resize(nBinsTotal,0);
        _mu.resize(nBinsTotal,0);
        _initialized.resize(nBinsTotal,false);
        _ds = separationBinning->getBinSize();
        _arcminToRad = 4*std::atan(1)/(60.*180.);
    }
    void addData(double value, double logLambda, double separation, double redshift) {
        // Lookup which (ll,sep,z) bin we are in.
        int logLambdaBin(_logLambdaBinning->getBinIndex(logLambda)),
            separationBin(_separationBinning->getBinIndex(separation)),
            redshiftBin(_redshiftBinning->getBinIndex(redshift));
        int index = (logLambdaBin*_nsep + separationBin)*_nz + redshiftBin;
        // Check that input (ll,sep,z) values correspond to bin centers.
        assert(std::fabs(logLambda-_logLambdaBinning->getBinCenter(logLambdaBin)) < 1e-6);
        assert(std::fabs(separation-_separationBinning->getBinCenter(separationBin)) < 1e-6);
        assert(std::fabs(redshift-_redshiftBinning->getBinCenter(redshiftBin)) < 1e-6);
        // Check that we have not already filled this bin.
        assert(!_initialized[index]);
        // Remember this bin.
        _data[index] = value;
        _initialized[index] = true;
        _index.push_back(index);
        _hasCov.push_back(false);
        // Calculate and save model observables for this bin.
        transform(logLambda,separation,redshift,_ds,_r3d[index],_mu[index]);
    }
    void transform(double ll, double sep, double z, double ds, double &r3d, double &mu) const {
        double ratio(std::exp(0.5*ll)),zp1(z+1);
        double z1(zp1/ratio-1), z2(zp1*ratio-1);
        double drLos = _cosmology->getLineOfSightComovingDistance(z2) -
            _cosmology->getLineOfSightComovingDistance(z1);
        // Calculate the geometrically weighted mean separation of this bin as
        // Integral[s^2,{s,smin,smax}]/Integral[s,{s,smin,smax}] = s + ds^2/(12*s)
        double swgt = sep + (ds*ds/12)/sep;
        double drPerp = _cosmology->getTransverseComovingScale(z)*(swgt*_arcminToRad);
        double rsq = drLos*drLos + drPerp*drPerp;
        r3d = std::sqrt(rsq);
        mu = std::abs(drLos)/r3d;
    }
    void addCovariance(int i, int j, double value) {
        assert(i >= 0 && i < getNData());
        // assert(j >= 0 && j < getNData());
        assert(i == j && value > 0);
        assert(_hasCov[i] == false);
        _cov[_index[i]] = value;
        _hasCov[i] = true;
    }
    int getSize() const { return _data.size(); }
    int getNData() const { return _index.size(); }
    int getNCov() const { return (int)std::count(_hasCov.begin(),_hasCov.end(),true); }
    int getIndex(int k) const { return _index[k]; }
    double getData(int index) const { return _data[index]; }
    double getVariance(int index) const { return _cov[index]; }
    double getRadius(int index) const { return _r3d[index]; }
    double getCosAngle(int index) const { return _mu[index]; }
    double getRedshift(int index) const { return _redshiftBinning->getBinCenter(index % _nz); }
    BinningPtr getLogLambdaBinning() const { return _logLambdaBinning; }
    BinningPtr getSeparationBinning() const { return _separationBinning; }
    BinningPtr getRedshiftBinning() const { return _redshiftBinning; }
private:
    BinningPtr _logLambdaBinning, _separationBinning, _redshiftBinning;
    cosmo::AbsHomogeneousUniversePtr _cosmology;
    std::vector<double> _data, _cov, _r3d, _mu;
    std::vector<bool> _initialized, _hasCov;
    std::vector<int> _index;
    int _ndata,_nsep,_nz;
    double _ds,_arcminToRad;
}; // LyaData

typedef boost::shared_ptr<LyaData> LyaDataPtr;

class LyaBaoModel {
public:
    LyaBaoModel(std::string const &fiducialName, std::string const &nowigglesName, double zref)
    : _zref(zref) {
        boost::format fileName("%s.%d.dat");
        _fid0 = load(boost::str(fileName % fiducialName % 0));
        _fid2 = load(boost::str(fileName % fiducialName % 2));
        _fid4 = load(boost::str(fileName % fiducialName % 4));
        _nw0 = load(boost::str(fileName % nowigglesName % 0));
        _nw2 = load(boost::str(fileName % nowigglesName % 2));
        _nw4 = load(boost::str(fileName % nowigglesName % 4));
        cosmo::CorrelationFunctionPtr fid0(new cosmo::CorrelationFunction(boost::bind(
            &lk::Interpolator::operator(),_fid0,_1)));
        cosmo::CorrelationFunctionPtr fid2(new cosmo::CorrelationFunction(boost::bind(
            &lk::Interpolator::operator(),_fid2,_1)));
        cosmo::CorrelationFunctionPtr fid4(new cosmo::CorrelationFunction(boost::bind(
            &lk::Interpolator::operator(),_fid4,_1)));
        cosmo::CorrelationFunctionPtr nw0(new cosmo::CorrelationFunction(boost::bind(
            &lk::Interpolator::operator(),_nw0,_1)));
        cosmo::CorrelationFunctionPtr nw2(new cosmo::CorrelationFunction(boost::bind(
            &lk::Interpolator::operator(),_nw2,_1)));
        cosmo::CorrelationFunctionPtr nw4(new cosmo::CorrelationFunction(boost::bind(
            &lk::Interpolator::operator(),_nw4,_1)));
        _fid.reset(new cosmo::RsdCorrelationFunction(fid0,fid2,fid4));
        _nw.reset(new cosmo::RsdCorrelationFunction(nw0,nw2,nw4));
    }
    double evaluate(double r, double mu, double z, lk::Parameters const &p) const {
        double alpha(p[0]), bias(p[1]), beta(p[2]), ampl(p[3]), scale(p[4]);
        double a1(p[5]), a2(p[6]), a3(p[7]);
        double zfactor = std::pow((1+z)/(1+_zref),alpha);
        _fid->setDistortion(beta);
        _nw->setDistortion(beta);
        double fid((*_fid)(r*scale,mu)), nw((*_nw)(r*scale,mu)); // scale cancels in mu
        double xi = ampl*(fid-nw)+nw;
        double broadband = 1e-1*a1/(r*r) + 1e-3*a2/r + 1e-5*a3;
        return bias*bias*zfactor*xi + broadband;
    }
private:
    lk::InterpolatorPtr load(std::string const &fileName) {
        std::vector<std::vector<double> > columns(2);
        std::ifstream in(fileName.c_str());
        lk::readVectors(in,columns);
        in.close();
        lk::InterpolatorPtr iptr(new lk::Interpolator(columns[0],columns[1],"cspline"));
        return iptr;
    }
    double _zref, _growth;
    lk::InterpolatorPtr _fid0, _fid2, _fid4, _nw0, _nw2, _nw4;
    boost::scoped_ptr<cosmo::RsdCorrelationFunction> _fid, _nw;
}; // LyaBaoModel

typedef boost::shared_ptr<LyaBaoModel> LyaBaoModelPtr;

typedef std::pair<double,double> ContourPoint;
typedef std::vector<ContourPoint> ContourPoints;

class Parameter {
public:
    Parameter(std::string const &name, double value, bool floating = false)
    : _name(name), _value(value), _floating(floating)
    { }
    void fix(double value) {
        _value = value;
        _floating = false;
    }
    void setValue(double value) { _value = value; }
    bool isFloating() const { return _floating; }
    double getValue() const { return _value; }
    std::string getName() const { return _name; }
private:
    std::string _name;
    double _value;
    bool _floating;
}; // Parameter

class LyaBaoLikelihood {
public:
    LyaBaoLikelihood(LyaDataPtr data, LyaBaoModelPtr model)
    : _data(data), _model(model),_errorScale(1) {
        assert(data);
        assert(model);
        _params.push_back(Parameter("Alpha",4.0,true));
        _params.push_back(Parameter("Bias",0.2,true));
        _params.push_back(Parameter("Beta",0.8,true));
        _params.push_back(Parameter("BAO Ampl",1,true));
        _params.push_back(Parameter("BAO Scale",1,true));
        _params.push_back(Parameter("BB a1",0,true));
        _params.push_back(Parameter("BB a2",0,true));
        _params.push_back(Parameter("BB a3",0,true));
    }
    void setErrorScale(double scale) {
        assert(scale > 0);
        _errorScale = scale;
    }
    double operator()(lk::Parameters const &params) {
        // Loop over the dataset bins.
        double nll(0);
        for(int k= 0; k < _data->getNData(); ++k) {
            int index(_data->getIndex(k));
            double r = _data->getRadius(index);
            double mu = _data->getCosAngle(index);
            double z = _data->getRedshift(index);
            double obs = _data->getData(index);
            double var = _data->getVariance(index);
            double pred = _model->evaluate(r,mu,z,params);
            // Update the chi2 = -log(L) for this bin
            double diff(obs-pred);
            nll += diff*diff/var;
        }
        // UP=0.5 is already hardcoded so we need a factor of 2 here since we are
        // calculating a chi-square. Apply an additional factor of _errorScale to
        // allow different error contours to be calculated.
        return 0.5*nll/_errorScale;
    }
    int getNPar() const { return _params.size(); }
    void initialize(lk::MinuitEngine::StatePtr initialState) {
        BOOST_FOREACH(Parameter const &param, _params) {
            double value(param.getValue());
            if(param.isFloating()) {
                initialState->Add(param.getName(),value,0.1*value); // error = 0.1*value
            }
            else {
                initialState->Add(param.getName(),value,0);
                initialState->Fix(param.getName());
            }
        }
    }
    void dump(std::string const &filename, lk::Parameters const &params,
    std::vector<ContourPoints> const &contourData,int oversampling = 10) {
        std::ofstream out(filename.c_str());
        // Dump binning info first
        BinningPtr llbins(_data->getLogLambdaBinning()), sepbins(_data->getSeparationBinning()),
            zbins(_data->getRedshiftBinning());
        out << llbins->getNBins() << ' ' << llbins->getLowEdge() << ' ' << llbins->getBinSize() << std::endl;
        out << sepbins->getNBins() << ' ' << sepbins->getLowEdge() << ' ' << sepbins->getBinSize() << std::endl;
        out << zbins->getNBins() << ' ' << zbins->getLowEdge() << ' ' << zbins->getBinSize() << std::endl;
        // Dump the number of data bins, the model oversampling factor, and the number of contour points.
        int ncontour = (0 == contourData.size()) ? 0 : contourData[0].size();
        out << _data->getNData() << ' ' << oversampling << ' ' << ncontour << std::endl;
        // Dump the number of parameters and their best-fit values.
        out << params.size();
        BOOST_FOREACH(double const &pValue, params) out << ' ' << pValue;
        out << std::endl;
        // Dump binned data and most recent pulls.
        for(int k= 0; k < _data->getNData(); ++k) {
            int index(_data->getIndex(k));
            double r = _data->getRadius(index);
            double mu = _data->getCosAngle(index);
            double z = _data->getRedshift(index);
            double obs = _data->getData(index);
            double var = _data->getVariance(index);
            double pred = _model->evaluate(r,mu,z,params);
            double pull = (obs-pred)/std::sqrt(var);
            out << index << ' ' << obs << ' ' << pull << std::endl;
        }
        // Dump oversampled model calculation.
        sepbins = oversampleBinning(*sepbins,oversampling);
        llbins = oversampleBinning(*llbins,oversampling);
        double r,mu,ds(sepbins->getBinSize());
        for(int iz = 0; iz < zbins->getNBins(); ++iz) {
            double z = zbins->getBinCenter(iz);
            for(int isep = 0; isep < sepbins->getNBins(); ++isep) {
                double sep = sepbins->getBinCenter(isep);
                for(int ill = 0; ill < llbins->getNBins(); ++ill) {
                    double ll = llbins->getBinCenter(ill);
                    _data->transform(ll,sep,z,ds,r,mu);
                    double pred = _model->evaluate(r,mu,z,params);
                    out << r << ' ' << pred << std::endl;
                }
            }
        }
        if(ncontour) {
            BOOST_FOREACH(ContourPoints const &points, contourData) {
                BOOST_FOREACH(ContourPoint const &point, points) {
                    out << point.first << ' ' << point.second << std::endl;
                }
            }
        }
        out.close();
    }
private:
    LyaDataPtr _data;
    LyaBaoModelPtr _model;
    std::vector<Parameter> _params;
    double _errorScale;
}; // LyaBaoLikelihood

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("BAO fitting");
    double OmegaLambda,OmegaMatter,zref,minll,dll,minsep,dsep,minz,dz;
    int nll,nsep,nz,ncontour;
    std::string fiducialName,nowigglesName,dataName,dumpName;
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.734),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0.266),
            "Present-day value of OmegaMatter or zero for 1-OmegaLambda.")
        ("fiducial", po::value<std::string>(&fiducialName)->default_value(""),
            "Fiducial correlation functions will be read from <name>.<ell>.dat with ell=0,2,4.")
        ("nowiggles", po::value<std::string>(&nowigglesName)->default_value(""),
            "No-wiggles correlation functions will be read from <name>.<ell>.dat with ell=0,2,4.")
        ("zref", po::value<double>(&zref)->default_value(2.25),
            "Reference redshift.")
        ("data", po::value<std::string>(&dataName)->default_value(""),
            "3D covariance data will be read from <data>.params and <data>.cov")
        ("minll", po::value<double>(&minll)->default_value(0.0002),
            "Minimum log(lam2/lam1).")
        ("dll", po::value<double>(&dll)->default_value(0.004),
            "log(lam2/lam1) binsize.")
        ("nll", po::value<int>(&nll)->default_value(14),
            "Maximum number of log(lam2/lam1) bins.")
        ("minsep", po::value<double>(&minsep)->default_value(0),
            "Minimum separation in arcmins.")
        ("dsep", po::value<double>(&dsep)->default_value(10),
            "Separation binsize in arcmins.")
        ("nsep", po::value<int>(&nsep)->default_value(14),
            "Maximum number of separation bins.")
        ("minz", po::value<double>(&minz)->default_value(1.7),
            "Minimum redshift.")
        ("dz", po::value<double>(&dz)->default_value(1.0),
            "Redshift binsize.")
        ("nz", po::value<int>(&nz)->default_value(2),
            "Maximum number of redshift bins.")
        ("dump", po::value<std::string>(&dumpName)->default_value(""),
            "Filename for dumping fit results.")
        ("ncontour",po::value<int>(&ncontour)->default_value(40),
            "Number of contour points to calculate in BAO parameters.")
        ("minos", "Runs MINOS to improve error estimates.")
        ;

    // Do the command line parsing now.
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose")), minos(vm.count("minos"));

    // Check for the required filename parameters.
    if(0 == dataName.length()) {
        std::cerr << "Missing required parameter --data." << std::endl;
        return -1;
    }
    if(0 == fiducialName.length()) {
        std::cerr << "Missing required parameter --fiducial." << std::endl;
        return -1;
    }
    if(0 == nowigglesName.length()) {
        std::cerr << "Missing required parameter --nowiggles." << std::endl;
        return -1;
    }

    // Initialize the cosmology calculations we will need.
    cosmo::AbsHomogeneousUniversePtr cosmology;
    LyaBaoModelPtr model;
    try {
        // Build the homogeneous cosmology we will use.
        if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
        cosmology.reset(new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));
        
         // Build our fit model from tabulated ell=0,2,4 correlation functions on disk.
         model.reset(new LyaBaoModel(fiducialName,nowigglesName,zref));

        if(verbose) std::cout << "Cosmology initialized." << std::endl;
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR during cosmology initialization:\n  " << e.what() << std::endl;
        return -2;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << "ERROR during cosmology initialization:\n  " << e.what() << std::endl;
        return -2;
    }
    
    // Load the data we will fit.
    LyaDataPtr data;
    try {
        // Initialize the (logLambda,separation,redshift) binning from command-line params.
        BinningPtr llBins(new Binning(nll,minll,dll)), sepBins(new Binning(nsep,minsep,dsep)),
            zBins(new Binning(nz,minz,dz));
        // Initialize the dataset we will fill.
        data.reset(new LyaData(llBins,sepBins,zBins,cosmology));
        // General stuff we will need for reading both files.
        std::string line;
        int lineNumber(0);
        // Capturing regexps for positive integer and signed floating-point constants.
        std::string ipat("(0|(?:[1-9][0-9]*))"),fpat("([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)");
        boost::match_results<std::string::const_iterator> what;
        // Loop over lines in the parameter file.
        std::string paramsName(dataName + ".params");
        std::ifstream paramsIn(paramsName.c_str());
        if(!paramsIn.good()) throw cosmo::RuntimeError("Unable to open " + paramsName);
        boost::regex paramPattern(
            boost::str(boost::format("\\s*%s\\s+%s\\s*\\| Lya covariance 3D \\(%s,%s,%s\\)\\s*")
            % fpat % fpat % fpat % fpat % fpat));
        while(paramsIn.good() && !paramsIn.eof()) {
            std::getline(paramsIn,line);
            if(paramsIn.eof()) break;
            if(!paramsIn.good()) {
                throw cosmo::RuntimeError("Unable to read line " + boost::lexical_cast<std::string>(lineNumber));
            }
            lineNumber++;
            // Parse this line with a regexp.
            if(!boost::regex_match(line,what,paramPattern)) {
                throw cosmo::RuntimeError("Badly formatted params line " +
                    boost::lexical_cast<std::string>(lineNumber) + ": '" + line + "'");
            }
            int nTokens(5);
            std::vector<double> token(nTokens);
            for(int tok = 0; tok < nTokens; ++tok) {
                token[tok] = boost::lexical_cast<double>(std::string(what[tok+1].first,what[tok+1].second));
            }
            // Add this bin to our dataset.
            if(0 != token[1]) throw cosmo::RuntimeError("Got unexpected non-zero token.");
            data->addData(token[0],token[2],token[3],token[4]);
        }
        paramsIn.close();
        if(verbose) {
            std::cout << "Read " << data->getNData() << " of " << data->getSize()
                << " data values from " << paramsName << std::endl;
        }
        // Loop over lines in the covariance file.
        std::string covName(dataName + ".cov");
        std::ifstream covIn(covName.c_str());
        if(!covIn.good()) throw cosmo::RuntimeError("Unable to open " + covName);
        boost::regex covPattern(boost::str(boost::format("\\s*%s\\s+%s\\s+%s\\s*") % ipat % ipat % fpat));
        lineNumber = 0;
        while(covIn.good() && !covIn.eof()) {
            std::getline(covIn,line);
            if(covIn.eof()) break;
            if(!covIn.good()) {
                throw cosmo::RuntimeError("Unable to read line " + boost::lexical_cast<std::string>(lineNumber));
            }
            lineNumber++;
            // Parse this line with a regexp.
            if(!boost::regex_match(line,what,covPattern)) {
                throw cosmo::RuntimeError("Badly formatted cov line " +
                    boost::lexical_cast<std::string>(lineNumber) + ": '" + line + "'");
            }
            int index1(boost::lexical_cast<int>(std::string(what[1].first,what[1].second)));
            int index2(boost::lexical_cast<int>(std::string(what[2].first,what[2].second)));
            double value(boost::lexical_cast<double>(std::string(what[3].first,what[3].second)));
            // Add this covariance to our dataset.
            data->addCovariance(index1,index2,value);
        }
        covIn.close();
        if(verbose) {
            std::cout << "Read " << data->getNCov() << " of " << data->getNData()
                << " diagonal covariance values from " << covName << std::endl;
        }
        assert(data->getNCov() == data->getNData());
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR while reading data:\n  " << e.what() << std::endl;
        return -2;
    }
    
    // Minimize the -log(Likelihood) function.
    try {
        lk::GradientCalculatorPtr gcptr;
        LyaBaoLikelihood nll(data,model);
        lk::FunctionPtr fptr(new lk::Function(boost::ref(nll)));

        int npar(nll.getNPar());
        lk::AbsEnginePtr engine = lk::getEngine("mn2::vmetric",fptr,gcptr,npar);
        lk::MinuitEngine &minuit = dynamic_cast<lk::MinuitEngine&>(*engine);        
        lk::MinuitEngine::StatePtr initialState(new ROOT::Minuit2::MnUserParameterState());
        nll.initialize(initialState);
        std::cout << *initialState;
        
        ROOT::Minuit2::MnStrategy strategy(1); // lo(0),med(1),hi(2)
        ROOT::Minuit2::MnMigrad fitter((ROOT::Minuit2::FCNBase const&)minuit,*initialState,strategy); 

        int maxfcn = 100*npar*npar;
        double edmtol = 0.1;
        ROOT::Minuit2::FunctionMinimum fmin = fitter(maxfcn,edmtol);
        
        if(minos) {
            ROOT::Minuit2::MnMinos minosError((ROOT::Minuit2::FCNBase const&)minuit,fmin,strategy);
            for(int ipar = 0; ipar < npar; ++ipar) {
                std::pair<double,double> error = minosError(ipar,maxfcn);
                std::cout << "MINOS error[" << ipar << "] = +" << error.second
                    << ' ' << error.first << std::endl;
            }
        }
        
        std::cout << fmin;
        std::cout << fmin.UserCovariance();
        std::cout << fmin.UserState().GlobalCC();
        
        std::vector<ContourPoints> contourData;
        if(ncontour > 0) {
            if(verbose) std::cout << "Calculating contours with " << ncontour << " points..." << std::endl;
            // 95% CL (see http://wwwasdoc.web.cern.ch/wwwasdoc/minuit/node33.html)
            // Calculate in mathematica using:
            // Solve[CDF[ChiSquareDistribution[2], x] == 0.95, x]
            nll.setErrorScale(5.99146);
            fmin = fitter(maxfcn,edmtol);
            ROOT::Minuit2::MnContours contours95((ROOT::Minuit2::FCNBase const&)minuit,fmin,strategy);
            // Parameter indices: 1=bias, 2=beta, 3=BAO amp, 4=BAO scale, 5=bband a1/10, 6=bband a2/1000
            contourData.push_back(contours95(5,6,ncontour));
            contourData.push_back(contours95(4,6,ncontour));
            contourData.push_back(contours95(1,6,ncontour));
            contourData.push_back(contours95(5,3,ncontour));
            contourData.push_back(contours95(4,3,ncontour));
            contourData.push_back(contours95(1,3,ncontour));            
            contourData.push_back(contours95(5,2,ncontour));
            contourData.push_back(contours95(4,2,ncontour));
            contourData.push_back(contours95(1,2,ncontour));
            // 68% CL
            nll.setErrorScale(2.29575);
            fmin = fitter(maxfcn,edmtol);
            ROOT::Minuit2::MnContours contours68((ROOT::Minuit2::FCNBase const&)minuit,fmin,strategy);
            contourData.push_back(contours68(5,6,ncontour));
            contourData.push_back(contours68(4,6,ncontour));
            contourData.push_back(contours68(1,6,ncontour));
            contourData.push_back(contours68(5,3,ncontour));
            contourData.push_back(contours68(4,3,ncontour));
            contourData.push_back(contours68(1,3,ncontour));            
            contourData.push_back(contours68(5,2,ncontour));
            contourData.push_back(contours68(4,2,ncontour));
            contourData.push_back(contours68(1,2,ncontour));
        }
        
        if(dumpName.length() > 0) {
            if(verbose) std::cout << "Dumping fit results to " << dumpName << std::endl;
            nll.dump(dumpName,fmin.UserParameters().Params(),contourData);
        }
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR during fit:\n  " << e.what() << std::endl;
        return -2;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << "ERROR during fit:\n  " << e.what() << std::endl;
        return -2;
    }

    // All done: normal exit.
    return 0;
}