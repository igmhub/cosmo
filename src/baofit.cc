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
#include <map>
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
        _ds2by12 = separationBinning->getBinSize()*separationBinning->getBinSize()/12;
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
        double ratio(std::exp(0.5*logLambda)),zp1(redshift+1);
        double z1(zp1/ratio-1), z2(zp1*ratio-1);
        double drLos = _cosmology->getLineOfSightComovingDistance(z2) -
            _cosmology->getLineOfSightComovingDistance(z1);
        // Calculate the geometrically weighted mean separation of this bin as
        // Integral[s^2,{s,smin,smax}]/Integral[s,{s,smin,smax}] = s + ds^2/(12*s)
        double swgt = separation + _ds2by12/separation;
        double drPerp = _cosmology->getTransverseComovingScale(redshift)*(swgt*_arcminToRad);
        double rsq = drLos*drLos + drPerp*drPerp;
        _r3d[index] = std::sqrt(rsq);
        _mu[index] = std::abs(drLos)/_r3d[index];
        // Agrees with Covariance3D::getRMuZ in ForestCovarianceParam.cpp
        /*
        std::cout << '(' << logLambda << ',' << separation << ',' << redshift << ") => ["
            << z1 << ',' << z2 << ',' << swgt << ';' << _drLos[index] << ','
            << _drPerp[index] << ',' << _mu[index] << "]\n";
        */
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
private:
    BinningPtr _logLambdaBinning, _separationBinning, _redshiftBinning;
    cosmo::AbsHomogeneousUniversePtr _cosmology;
    std::vector<double> _data, _cov, _r3d, _mu;
    std::vector<bool> _initialized, _hasCov;
    std::vector<int> _index;
    int _ndata,_nsep,_nz;
    double _ds2by12,_arcminToRad;
}; // LyaData

typedef boost::shared_ptr<LyaData> LyaDataPtr;

class Parameter {
public:
    Parameter(double value, bool floating = false)
    : _value(value), _floating(floating)
    { }
    void fix(double value) {
        _value = value;
        _floating = false;
    }
    void setValue(double value) { _value = value; }
    bool isFloating() const { return _floating; }
    double getValue() const { return _value; }
private:
    double _value;
    bool _floating;
}; // Parameter

typedef std::pair<std::string,Parameter> NamedParameter;
typedef std::map<std::string,Parameter> ParameterMap;

class LyaBaoLikelihood {
public:
    LyaBaoLikelihood(LyaDataPtr data, BaoFitPowerPtr power, double zref, double growth,
    double rmin, double rmax, int nr) : _data(data), _power(power), _zref(zref), _growth(growth),
    _pptr(new cosmo::PowerSpectrum(boost::ref(*_power))), _xi(_pptr,rmin,rmax,nr),
    _rmin(rmin), _rmax(rmax)
    {
        assert(data);
        assert(power);
        _params.insert(NamedParameter("Alpha",Parameter(4.0,true)));
        _params.insert(NamedParameter("Bias",Parameter(0.2,true)));
        _params.insert(NamedParameter("Beta",Parameter(0.8,true)));
        _params.insert(NamedParameter("BAO Ampl",Parameter(1,false)));
        _params.insert(NamedParameter("BAO Scale",Parameter(1,false)));
        _params.insert(NamedParameter("BAO Sigma",Parameter(0,false))); // in Mpc/h
    }
    double operator()(lk::Parameters const &params) {
        // Update the values of each parameter.
        lk::Parameters::const_iterator nextValue(params.begin());
        ParameterMap::iterator iter;
        for(iter = _params.begin(); iter != _params.end(); ++iter) {
            iter->second.setValue(*nextValue++);
        }
        // Transfer parameter values to the appropriate models.
        double alpha(_params.find("Alpha")->second.getValue());
        double bias(_params.find("Bias")->second.getValue());
        _xi.setDistortion(_params.find("Beta")->second.getValue());
        _power->setAmplitude(_params.find("BAO Ampl")->second.getValue());
        _power->setScale(_params.find("BAO Scale")->second.getValue());
        _power->setSigma(_params.find("BAO Sigma")->second.getValue());
        // Loop over the dataset bins.
        double nll(0);
        double biasSq(bias*bias);
        for(int k= 0; k < _data->getNData(); ++k) {
            int index(_data->getIndex(k));
            double r = _data->getRadius(index);
            if(r < _rmin || r >= _rmax) continue;
            double mu = _data->getCosAngle(index);
            double z = _data->getRedshift(index);
            double zfactor = _growth*std::pow((1+z)/(1+_zref),alpha);
            double pred = biasSq*zfactor*_xi(r,mu);
            double obs = _data->getData(index);
            double var = _data->getVariance(index);
            // Update the chi2 = -log(L) for this bin
            double diff(obs-pred);
            nll += diff*diff/var;
            /**
            std::cout << index << ' ' << r << ' ' << mu << ' ' << z << " => "
                << pred << ' ' << obs << ' ' << err << std::endl;
            **/
        }
        return 0.5*nll; // convert chi2 into -log(L) to match UP=1
    }
    int getNPar() const { return _params.size(); }
    void initialize(lk::MinuitEngine::StatePtr initialState) {
        BOOST_FOREACH(NamedParameter const &np, _params) {
            std::string const &name(np.first);
            double value(np.second.getValue());
            if(np.second.isFloating()) {
                initialState->Add(name,value,0.1*value); // error = 0.1*value
            }
            else {
                initialState->Add(name,value,0);
                initialState->Fix(name);
            }
        }
    }
    lk::Parameters getInitialValues() const {
        lk::Parameters initial;
        BOOST_FOREACH(NamedParameter const &np, _params) {
            if(np.second.isFloating()) initial.push_back(np.second.getValue());
        }
        return initial;
    }
    lk::Parameters getInitialErrors() const {
        lk::Parameters errors;
        BOOST_FOREACH(NamedParameter const &np, _params) {
            if(np.second.isFloating()) errors.push_back(0.1*np.second.getValue());
        }
        return errors;
    }
private:
    LyaDataPtr _data;
    BaoFitPowerPtr _power;
    cosmo::PowerSpectrumPtr _pptr;
    cosmo::RsdPowerSpectrumCorrelationFunction _xi;
    double _zref, _growth, _rmin, _rmax;
    ParameterMap _params;
}; // LyaBaoLikelihood

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("BAO fitting");
    double OmegaLambda,OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,spectralIndex,sigma8,zref;
    double minll,dll,minsep,dsep,minz,dz,rmin,rmax;
    int nll,nsep,nz,nr;
    std::string dataName;
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.734),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0.266),
            "Present-day value of OmegaMatter or zero for 1-OmegaLambda.")
        ("omega-baryon", po::value<double>(&OmegaBaryon)->default_value(0.0449),
            "Present-day value of OmegaBaryon, must be <= OmegaMatter.")
        ("hubble-constant", po::value<double>(&hubbleConstant)->default_value(0.710),
            "Present-day value of the Hubble parameter h = H0/(100 km/s/Mpc).")
        ("cmb-temp", po::value<double>(&cmbTemp)->default_value(2.725),
            "Present-day temperature of the cosmic microwave background in Kelvin.")
        ("spectral-index", po::value<double>(&spectralIndex)->default_value(1),
            "Power exponent of primordial fluctuations.")
        ("sigma8", po::value<double>(&sigma8)->default_value(0.801),
            "Power will be normalized to this value.")
        ("zref", po::value<double>(&zref)->default_value(2.25),
            "Reference redshift.")
        ("data", po::value<std::string>(&dataName)->default_value(""),
            "3D covariance data will be read from <data>.params and <data>.cov")
        ("minll", po::value<double>(&minll)->default_value(0.0002),
            "Minimum log(lam2/lam1).")
        ("dll", po::value<double>(&dll)->default_value(0.004),
            "log(lam2/lam1) binsize.")
        ("nll", po::value<int>(&nll)->default_value(99),
            "Maximum number of log(lam2/lam1) bins.")
        ("minsep", po::value<double>(&minsep)->default_value(0),
            "Minimum separation in arcmins.")
        ("dsep", po::value<double>(&dsep)->default_value(10),
            "Separation binsize in arcmins.")
        ("nsep", po::value<int>(&nsep)->default_value(18),
            "Maximum number of separation bins.")
        ("minz", po::value<double>(&minz)->default_value(1.7),
            "Minimum redshift.")
        ("dz", po::value<double>(&dz)->default_value(1.0),
            "Redshift binsize.")
        ("nz", po::value<int>(&nz)->default_value(2),
            "Maximum number of redshift bins.")
        ("rmin", po::value<double>(&rmin)->default_value(10),
            "Minimum value of 3D separation of fit in Mpc/h.")
        ("rmax", po::value<double>(&rmax)->default_value(170),
            "Maximum value of 3D separation of fit in Mpc/h.")
        ("nr", po::value<int>(&nr)->default_value(1024),
            "Number of interpolation points to use in 3D separation.")
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
    bool verbose(vm.count("verbose"));
    // Check for the required data name.
    if(0 == dataName.length()) {
        std::cerr << "Missing required parameter --data." << std::endl;
        return -1;
    }

    // Initialize the cosmology calculations we will need.
    cosmo::AbsHomogeneousUniversePtr cosmology;
    BaoFitPowerPtr power;
    double redshiftGrowthFactor(1);
    try {
        // Build the homogeneous cosmology we will use.
        if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
        cosmology.reset(new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));
        
        // Calculate the redshift growth factor to apply for the reference redshift.
         double ratio = cosmology->getGrowthFunction(zref)/cosmology->getGrowthFunction(0);
         redshiftGrowthFactor = ratio*ratio;
         if(verbose) {
             std::cout << "Redshift growth factor for zref = " << zref << " is " << redshiftGrowthFactor << std::endl;
         }
    
        // Build fiducial and "no-wiggles" Eisenstein & Hu models.
        boost::shared_ptr<cosmo::BaryonPerturbations>
            baryons(new cosmo::BaryonPerturbations(
                OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,cosmo::BaryonPerturbations::ShiftedOscillation)),
            nowiggles(new cosmo::BaryonPerturbations(
                OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,cosmo::BaryonPerturbations::NoOscillation));
        
        // Dump some info about the fiducial model if requested.
        if(verbose) {
            std::cout << "z(eq) = " << baryons->getMatterRadiationEqualityRedshift() << std::endl;
            std::cout << "k(eq) = " << baryons->getMatterRadiationEqualityScale() << " /(Mpc/h)"
                << std::endl;
            std::cout << "sound horizon = " << baryons->getSoundHorizon() << " Mpc/h at z(drag) = "
                << baryons->getDragEpoch() << std::endl;
            std::cout << "Silk damping scale = " << baryons->getSilkDampingScale() << " /(Mpc/h)"
                << std::endl;
        }

        // Make shareable pointers to the matter transfer functions of these models.
        cosmo::TransferFunctionPtr
            baryonsTransferPtr(new cosmo::TransferFunction(boost::bind(
                &cosmo::BaryonPerturbations::getMatterTransfer,baryons,_1))),
            nowigglesTransferPtr(new cosmo::TransferFunction(boost::bind(
                &cosmo::BaryonPerturbations::getMatterTransfer,nowiggles,_1)));    

        // Build the corresponding power spectra.
        boost::shared_ptr<cosmo::TransferFunctionPowerSpectrum>
            baryonsPower(new cosmo::TransferFunctionPowerSpectrum(baryonsTransferPtr,spectralIndex)),
            nowigglesPower(new cosmo::TransferFunctionPowerSpectrum(nowigglesTransferPtr,spectralIndex));
        // Normalize the fiducial model to sigma8, and use the same value of deltaH for the nowiggles model.
        baryonsPower->setSigma(sigma8,8);
        nowigglesPower->setDeltaH(baryonsPower->getDeltaH());

        // Make shareable power spectrum pointers.
        cosmo::PowerSpectrumPtr
            baryonsPowerPtr(new cosmo::PowerSpectrum(boost::bind(
                &cosmo::TransferFunctionPowerSpectrum::operator(),baryonsPower,_1))),
            nowigglesPowerPtr(new cosmo::PowerSpectrum(boost::bind(
                &cosmo::TransferFunctionPowerSpectrum::operator(),nowigglesPower,_1)));

        // Build a hybrid power spectrum that combines the fiducial and nowiggles models.
        power.reset(new BaoFitPower(baryonsPowerPtr,nowigglesPowerPtr));
        
        if(verbose) std::cout << "Cosmology initialized." << std::endl;
    }
    catch(cosmo::RuntimeError const &e) {
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
        LyaBaoLikelihood nll(data,power,zref,redshiftGrowthFactor,rmin,rmax,nr);
        lk::FunctionPtr fptr(new lk::Function(nll));

        int npar(nll.getNPar());
        lk::AbsEnginePtr engine = lk::getEngine("mn2::vmetric",fptr,gcptr,npar);
        lk::MinuitEngine &minuit = dynamic_cast<lk::MinuitEngine&>(*engine);        
        lk::MinuitEngine::StatePtr initialState(new ROOT::Minuit2::MnUserParameterState());
        nll.initialize(initialState);
        std::cout << *initialState;
        
        ROOT::Minuit2::MnMigrad fitter((ROOT::Minuit2::FCNBase const&)(minuit),*initialState,
            ROOT::Minuit2::MnStrategy(1)); // lo(0),med(1),hi(2)

        int maxfcn = 100*npar*npar;
        double edmtol = 0.1;
        ROOT::Minuit2::FunctionMinimum min = fitter(maxfcn,edmtol);
        std::cout << min;
        std::cout << min.UserCovariance();
        std::cout << min.UserState().GlobalCC();

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