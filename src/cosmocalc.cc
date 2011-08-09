// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/cosmo.h"

#include "boost/program_options.hpp"

#include <iostream>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology calculator");
    double OmegaLambda,OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,zval,kval;
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.728),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0),
            "Present-day value of OmegaMatter or zero for 1-OmegaLambda.")
        ("omega-baryon", po::value<double>(&OmegaBaryon)->default_value(0.0456),
            "Present-day value of OmegaBaryon, must be <= OmegaMatter.")
        ("hubble-constant", po::value<double>(&hubbleConstant)->default_value(0.704),
            "Present-day value of the Hubble parameter h = H0/(100 km/s/Mpc).")
        ("cmb-temp", po::value<double>(&cmbTemp)->default_value(2.725),
            "Present-day temperature of the cosmic microwave background in Kelvin.")
        ("redshift,z", po::value<double>(&zval)->default_value(1),
            "Emitter redshift.")
        ("wavenumber,k", po::value<double>(&kval)->default_value(0.1),
            "Perturbation wavenumber in 1/(Mpc/h).")
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

    if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
    cosmo::AbsHomogeneousUniversePtr cosmology(
        new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));
    
    std::cout << "z = " << zval << std::endl;
    std::cout << "D(z) = " << cosmology->getLineOfSightComovingDistance(zval) << " Mpc/h"
        << std::endl;
    std::cout << "D1(z) = " << 2.5*OmegaMatter*cosmology->getGrowthFunction(zval)
        << std::endl;
    
    cosmo::BaryonPerturbations baryons(OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp);

    std::cout << "z(eq) = " << baryons.getMatterRadiationEqualityRedshift() << std::endl;
    std::cout << "k(eq) = " << baryons.getMatterRadiationEqualityScale() << " /(Mpc/h)"
        << std::endl;
    std::cout << "sound horizon = " << baryons.getSoundHorizon() << " Mpc/h at z(drag) = "
        << baryons.getDragEpoch() << std::endl;
    std::cout << "Silk damping scale = " << baryons.getSilkDampingScale() << " /(Mpc/h)"
        << std::endl;

    double Tfc,Tfb,Tf;
    baryons.calculateTransferFunctions(kval,Tfc,Tfb,Tf);
    std::cout << "k = " << kval << " /(Mpc/h)" << std::endl;
    std::cout << "Tf(cmb,k) = " << Tfc << std::endl;
    std::cout << "Tf(baryon,k) = " << Tfb << std::endl;
    std::cout << "Tf(full,k) = " << Tf << std::endl;

    return 0;
}