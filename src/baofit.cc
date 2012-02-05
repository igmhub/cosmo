// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/cosmo.h"

#include "boost/program_options.hpp"
#include "boost/bind.hpp"
#include "boost/ref.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "boost/format.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

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
    LyaData(BinningPtr logLambdaBinning, BinningPtr separationBinning, BinningPtr redshiftBinning) :
    _logLambdaBinning(logLambdaBinning), _separationBinning(separationBinning), _redshiftBinning(redshiftBinning)
    {
        assert(logLambdaBinning);
        assert(separationBinning);
        assert(redshiftBinning);
        _nsep = separationBinning->getNBins();
        _nz = redshiftBinning->getNBins();
        int nBinsTotal = logLambdaBinning->getNBins()*_nsep*_nz;
        _data.resize(nBinsTotal,0);
        _cov.resize(nBinsTotal,0);
        _initialized.resize(nBinsTotal,false);
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
    }
    void addCovariance(int i, int j, double value) {
        assert(i >= 0 && i < getNData());
        // assert(j >= 0 && j < getNData());
        assert(i == j);
        assert(_hasCov[i] == false);
        _cov[_index[i]] = value;
        _hasCov[i] = true;
    }
    int getSize() const { return _data.size(); }
    int getNData() const { return _index.size(); }
    int getNCov() const { return (int)std::count(_hasCov.begin(),_hasCov.end(),true); }
private:
    BinningPtr _logLambdaBinning, _separationBinning, _redshiftBinning;
    std::vector<double> _data, _cov;
    std::vector<bool> _initialized, _hasCov;
    std::vector<int> _index;
    int _ndata,_nsep,_nz;
}; // LyaData

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("BAO fitting");
    double OmegaLambda,OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,spectralIndex,sigma8,zref;
    double minll,dll,minsep,dsep,minz,dz;
    int nll,nsep,nz;
    std::string dataName;
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.734),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0.226),
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
    try {
        // Build the homogeneous cosmology we will use.
        if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
        cosmo::AbsHomogeneousUniversePtr cosmology(new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));
    
        // Build fiducial and "no-wiggles" Eisenstein & Hu models.
        cosmo::BaryonPerturbations
            baryons(OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,cosmo::BaryonPerturbations::ShiftedOscillation),
            nowiggles(OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,cosmo::BaryonPerturbations::NoOscillation);
        
        // Dump some info about the fiducial model if requested.
        if(verbose) {
            std::cout << "z(eq) = " << baryons.getMatterRadiationEqualityRedshift() << std::endl;
            std::cout << "k(eq) = " << baryons.getMatterRadiationEqualityScale() << " /(Mpc/h)"
                << std::endl;
            std::cout << "sound horizon = " << baryons.getSoundHorizon() << " Mpc/h at z(drag) = "
                << baryons.getDragEpoch() << std::endl;
            std::cout << "Silk damping scale = " << baryons.getSilkDampingScale() << " /(Mpc/h)"
                << std::endl;
        }

        // Make shareable pointers to the matter transfer functions of these models.
        cosmo::TransferFunctionPtr
            baryonsTransferPtr(new cosmo::TransferFunction(boost::bind(
                &cosmo::BaryonPerturbations::getMatterTransfer,&baryons,_1))),
            nowigglesTransferPtr(new cosmo::TransferFunction(boost::bind(
                &cosmo::BaryonPerturbations::getMatterTransfer,&nowiggles,_1)));    

        // Build the corresponding power spectra.
        cosmo::TransferFunctionPowerSpectrum
            baryonsPower(baryonsTransferPtr,spectralIndex),
            nowigglesPower(nowigglesTransferPtr,spectralIndex);
        // Normalize the fiducial model to sigma8, and use the same value of deltaH for the nowiggles model.
        baryonsPower.setSigma(sigma8,8);
        nowigglesPower.setDeltaH(baryonsPower.getDeltaH());

        // Make shareable power spectrum pointers.
        cosmo::PowerSpectrumPtr
            baryonsPowerPtr(new cosmo::PowerSpectrum(boost::ref(baryonsPower))),
            nowigglesPowerPtr(new cosmo::PowerSpectrum(boost::ref(nowigglesPower)));

        // Build a hybrid power spectrum that combines the fiducial and nowiggles models.
        BaoFitPower hybridPower(baryonsPowerPtr,nowigglesPowerPtr);
        cosmo::PowerSpectrumPtr hybridPowerPtr(new cosmo::PowerSpectrum(boost::ref(hybridPower)));
        
        if(verbose) std::cout << "Cosmology initialized." << std::endl;
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR during cosmology initialization:\n  " << e.what() << std::endl;
        return -2;
    }
    
    // Load the data we will fit.
    try {
        // Initialize the (logLambda,separation,redshift) binning from command-line params.
        BinningPtr llBins(new Binning(nll,minll,dll)), sepBins(new Binning(nsep,minsep,dsep)),
            zBins(new Binning(nz,minz,dz));
        // Initialize the dataset we will fill.
        LyaData dataset(llBins,sepBins,zBins);
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
            dataset.addData(token[0],token[2],token[3],token[4]);
        }
        paramsIn.close();
        if(verbose) {
            std::cout << "Read " << dataset.getNData() << " of " << dataset.getSize()
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
            dataset.addCovariance(index1,index2,value);
        }
        covIn.close();
        if(verbose) {
            std::cout << "Read " << dataset.getNCov() << " of " << dataset.getNData()
                << " diagonal covariance values from " << covName << std::endl;
        }
        assert(dataset.getNCov() == dataset.getNData());
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR while reading data:\n  " << e.what() << std::endl;
        return -2;
    }

    // All done: normal exit.
    return 0;
}