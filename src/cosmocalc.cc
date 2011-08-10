// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

// Reproduce bottom-left plot of Fig.3 in astro-ph/9709112 using:
// cosmocalc --omega-matter 0.2 --omega-baryon 0.1 --hubble-constant 0.5 --cmb-temp 2.728 \
//   --kmin 0.001 --kmax 1 --nk 500 --save-transfer fig3.dat

#include "cosmo/cosmo.h"

#include "boost/program_options.hpp"
#include "boost/bind.hpp"
#include "boost/ref.hpp"

#include <fstream>
#include <iostream>
#include <cmath>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology calculator");
    double OmegaLambda,OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,spectralIndex,
        zval,kval,kmin,kmax;
    int nk;
    std::string saveTransferFile;
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
        ("spectral-index", po::value<double>(&spectralIndex)->default_value(1),
            "Power exponent of primordial fluctuations.")
        ("redshift,z", po::value<double>(&zval)->default_value(1),
            "Emitter redshift.")
        ("wavenumber,k", po::value<double>(&kval)->default_value(0.1),
            "Perturbation wavenumber in 1/(Mpc/h).")
        ("save-transfer", po::value<std::string>(&saveTransferFile)->default_value(""),
            "Saves the matter transfer function to the specified filename.")
        ("kmin", po::value<double>(&kmin)->default_value(0.001),
            "Minimum wavenumber in 1/(Mpc/h) for tabulating transfer function.")
        ("kmax", po::value<double>(&kmax)->default_value(100.),
            "Maximum wavenumber in 1/(Mpc/h) for tabulating transfer function.")
        ("nk", po::value<int>(&nk)->default_value(100),
            "Number of logarithmic steps to use for tabulating transfer function.")
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
    std::cout << "DM(z) = " << cosmology->getTransverseComovingScale(zval) << " Mpc/h/rad"
        << std::endl;
    double tL(cosmology->getLookbackTime(zval));
    double conv(1e9*86400*365.25);
    std::cout << "t(lookback,z) = " << tL << " secs/h = " << tL/conv*hubbleConstant
        << " Gyr" << std::endl;
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
    
    if(0 < saveTransferFile.length()) {
        std::ofstream out(saveTransferFile.c_str());
        double kratio(std::pow(kmax/kmin,1/(nk-1.)));
        for(int i = 0; i < nk; ++i) {
            double k(kmin*std::pow(kratio,i));
            double Tf(baryons.getMatterTransfer(k));
            out << k << ' ' << Tf << std::endl;
        }
        out.close();
    }

    cosmo::TransferFunctionPtr transfer(new cosmo::TransferFunction(boost::bind(
        &cosmo::BaryonPerturbations::getMatterTransfer,&baryons,_1)));

    // Use COBE normalization for n=1
    double deltaH(1.94e-5*std::pow(OmegaMatter,-0.785-0.05*std::log(OmegaMatter)));
    std::cout << "deltaH = " << deltaH << std::endl;

    cosmo::TransferFunctionPowerSpectrum transferPower(transfer,spectralIndex,deltaH);
    cosmo::PowerSpectrumPtr power(new cosmo::PowerSpectrum(boost::ref(transferPower)));
    double sig8pred(0.5*std::pow(OmegaMatter,-0.65));
    
    std::cout << "sigma(8 Mpc/h) = " << cosmo::getRmsAmplitude(power,8)
        << " (pred = " << sig8pred << ")" << std::endl;
        
    // Calculate the Gaussian RMS amplitude on the Jean's length appropriate for
    // QSO spectra, evolved for z = 3.
    double rQSO(0.0416/std::sqrt(OmegaMatter)); // in Mpc/h
    double evol(cosmology->getGrowthFunction(3)/cosmology->getGrowthFunction(0));
    double sigmaQSO(cosmo::getRmsAmplitude(power,rQSO,true));

    std::cout << "rQSO = " << rQSO << " Mpc/h, sigmaQSO(z=0) = " << sigmaQSO
        << ", sigmaQSO(z=3) = " << sigmaQSO*evol << std::endl;

    return 0;
}