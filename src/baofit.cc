// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/cosmo.h"

#include "boost/program_options.hpp"
#include "boost/bind.hpp"
#include "boost/ref.hpp"

#include <fstream>
#include <iostream>
#include <cmath>

namespace po = boost::program_options;

class BaoFitPower {
public:
    BaoFitPower(double amplitude, double scale, double sigma,
        cosmo::PowerSpectrumPtr full, cosmo::PowerSpectrumPtr nowiggles)
    : _amplitude(amplitude), _scale(scale), _scale4(scale*scale*scale*scale), _sigsq(sigma*sigma),
        _full(full), _nowiggles(nowiggles)
    { }
    double operator()(double k) const {
        double ak(k/_scale), smooth(std::exp(-ak*ak*_sigsq/2));
        double fullPower = (*_full)(ak), nowigglesPower = (*_nowiggles)(ak);
        return _scale4*(_amplitude*smooth*(fullPower - nowigglesPower) + nowigglesPower);
    }
private:
    double _amplitude, _scale, _scale4, _sigsq;
    cosmo::PowerSpectrumPtr _full, _nowiggles;
}; // BaoFitPower

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology calculator");
    double OmegaLambda,OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,spectralIndex,sigma8,zref;
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
        ;

    // do the command line parsing now
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

    // Build the homogeneous cosmology we will use.
    if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
    cosmo::AbsHomogeneousUniversePtr cosmology(new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));
    
    // Build fiducial and "no-wiggles" Eisenstein & Hu models.
    cosmo::BaryonPerturbations
        baryons(OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,cosmo::BaryonPerturbations::ShiftedOscillation),
        nowiggles(OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,cosmo::BaryonPerturbations::NoOscillation);

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

    //!!cosmo::PowerSpectrumPtr power(new cosmo::PowerSpectrum(boost::ref(transferPower)));

    std::cout << "z(eq) = " << baryons.getMatterRadiationEqualityRedshift() << std::endl;
    std::cout << "k(eq) = " << baryons.getMatterRadiationEqualityScale() << " /(Mpc/h)"
        << std::endl;
    std::cout << "sound horizon = " << baryons.getSoundHorizon() << " Mpc/h at z(drag) = "
        << baryons.getDragEpoch() << std::endl;
    std::cout << "Silk damping scale = " << baryons.getSilkDampingScale() << " /(Mpc/h)"
        << std::endl;

    return 0;
}