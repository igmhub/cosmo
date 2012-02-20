// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

// Reproduce bottom-left plot of Fig.3 in astro-ph/9709112 using:
// cosmocalc --omega-matter 0.2 --omega-baryon 0.1 --hubble-constant 0.5 --cmb-temp 2.728 \
//   --kmin 0.001 --kmax 1 --nk 500 --save-power fig3.dat

// Reproduce Fig.1 of JMLG paper draft (needs an extra factor of pi/2 ??)
// ./cosmocalc --omega-baryon 0.044 --omega-matter 0.27 --omega-lambda 0.73 \
//   --hubble-constant 0.71 --save-power xfer.dat -r 0.1 --kmax 1

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/bind.hpp"

#include <fstream>
#include <iostream>
#include <cmath>

namespace po = boost::program_options;
namespace lk = likely;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology calculator");
    double OmegaLambda,OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,spectralIndex,sigma8,
        zval,kval,kmin,kmax,r1d,rmin,rmax,baoAmplitude,baoSigma,baoScale;
    double bbandP,bbandCoef,bbandKmin,bbandRmin,bbandR0,bbandVar;
    int nk,nr;
    std::string loadPowerFile,savePowerFile,saveCorrelationFile;
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
        ("sigma8", po::value<double>(&sigma8)->default_value(0),
            "Power will be normalized to this value (default is COBE normalization).")
        ("redshift,z", po::value<double>(&zval)->default_value(1),
            "Emitter redshift.")
        ("wavenumber,k", po::value<double>(&kval)->default_value(0.1),
            "Wavenumber in h/Mpc to use for verbose output.")
        ("load-power", po::value<std::string>(&loadPowerFile)->default_value(""),
            "Reads k,P(k) values (in h/Mpc units) to interpolate from the specified filename.")
        ("save-power", po::value<std::string>(&savePowerFile)->default_value(""),
            "Saves the matter power spectrum to the specified filename.")
        ("power1d", "Adds 1D power spectrum to save-power output file.")
        ("r1d,r", po::value<double>(&r1d)->default_value(0.04),
            "Radius for calculating 1D power spectrum in Mpc/h.")
        ("kmin", po::value<double>(&kmin)->default_value(0.001),
            "Minimum wavenumber in 1/(Mpc/h) for tabulating transfer function.")
        ("kmax", po::value<double>(&kmax)->default_value(100.),
            "Maximum wavenumber in 1/(Mpc/h) for tabulating transfer function.")
        ("nk", po::value<int>(&nk)->default_value(100),
            "Number of logarithmic steps to use for tabulating transfer function.")
        ("save-correlation", po::value<std::string>(&saveCorrelationFile)->default_value(""),
            "Saves the matter correlation function to the specified filename.")
        ("rmin", po::value<double>(&rmin)->default_value(0.01),
            "Minimum radius in (Mpc/h) for tabulating correlation function.")
        ("rmax", po::value<double>(&rmax)->default_value(1000.),
            "Maximum radius in (Mpc/h) for tabulating correlation function.")
        ("nr", po::value<int>(&nr)->default_value(100),
            "Number of logarithmic steps to use for tabulating correlation function.")
        ("rlog", "Use log spaced r-values for saved correlation function (default is linear).")
        ("quad", "Calculates the quadrupole (l=2) correlation function (default is monopole).")
        ("hexa", "Calculates the hexedacapole (l=4) correlation function (default is monopole).")
        ("no-wiggles", "Calculates the power spectrum without baryon acoustic oscillations.")
        ("periodic-wiggles", "Calculates the power spectrum with periodic acoustic oscillations.")
        ("bao-amplitude", po::value<double>(&baoAmplitude)->default_value(1),
            "Amplitude of baryon acoustic oscillations relative to fiducial model.")
        ("bao-sigma", po::value<double>(&baoSigma)->default_value(0),
            "Gaussian smearing of BAO correlation function peak in Mpc/h relative to fiducial model.")
        ("bao-scale", po::value<double>(&baoScale)->default_value(1),
            "Rescaling of wavenumber relative to fiducial model (>1 means larger acoustic scale)")
        ("broadband-p", po::value<double>(&bbandP)->default_value(0),
            "Broadband exponent (1/k)^p.")
        ("broadband-coef", po::value<double>(&bbandCoef)->default_value(0),
            "Coefficient of broadband power.")
        ("broadband-kmin", po::value<double>(&bbandKmin)->default_value(0),
            "Small k cutoff for regulating broadband power in h/Mpc.")
        ("broadband-rmin", po::value<double>(&bbandRmin)->default_value(0),
            "Small r (large k) cutoff for regulating broadband power in Mpc/h.")
        ("broadband-r0", po::value<double>(&bbandR0)->default_value(8),
            "Scale for fixing 'natural normalization' of broadband power.")
        ("broadband-var", po::value<double>(&bbandVar)->default_value(0.1),
            "Variance within r0 for fixing 'natural normalization' of broadband power.")
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
    bool verbose(vm.count("verbose")), power1d(vm.count("power1d")), rlog(vm.count("rlog")),
        quad(vm.count("quad")), hexa(vm.count("hexa")), noWiggles(vm.count("no-wiggles")),
        periodicWiggles(vm.count("periodic-wiggles"));

    // Process the multipole flags.
    if(quad && hexa) {
        std::cerr << "Cannot request both quad (l=2) and hexa (l=4) for correlation function output."
            << std::endl;
        return -1;
    }
    cosmo::PowerSpectrumCorrelationFunction::Multipole
        multipole(cosmo::PowerSpectrumCorrelationFunction::Monopole);
    if(quad) multipole = cosmo::PowerSpectrumCorrelationFunction::Quadrupole;
    if(hexa) multipole = cosmo::PowerSpectrumCorrelationFunction::Hexadecapole;

    // Process the wiggle flags.
    if(vm.count("no-wiggles")+vm.count("periodic-wiggles")+vm.count("bao-fit") > 1) {
        std::cerr << "Specify at most one of no-wiggles, periodic-wiggles, bao-fit options."
            << std::endl;
        return -1;
    }
    cosmo::BaryonPerturbations::BaoOption baoOption(cosmo::BaryonPerturbations::ShiftedOscillation);
    if(noWiggles) baoOption = cosmo::BaryonPerturbations::NoOscillation;
    if(periodicWiggles) baoOption = cosmo::BaryonPerturbations::PeriodicOscillation;

    // Build the homogeneous cosmology we will use.
    if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
    cosmo::AbsHomogeneousUniversePtr cosmology(
        new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));
    double growthFactor = cosmology->getGrowthFunction(zval)/cosmology->getGrowthFunction(0);
    if(verbose) {
        // Print homogeneous cosmology info.
        std::cout << "curvature = " << cosmology->getCurvature() << std::endl;    
        std::cout << "At z = " << zval << ':' << std::endl;
        std::cout << "  H(z)/H0 = " << cosmology->getHubbleFunction(zval) << std::endl;
        std::cout << "  Radial D(z) = " << cosmology->getLineOfSightComovingDistance(zval)
            << " Mpc/h" << std::endl;
        std::cout << "  Transverse DA(z) = " << cosmology->getTransverseComovingScale(zval)
            << " Mpc/h/rad" << std::endl;
        double tL(cosmology->getLookbackTime(zval));
        double conv(1e9*86400*365.25);
        std::cout << "  t(lookback,z) = " << tL << " secs/h = " << tL/conv*hubbleConstant
            << " Gyr/h" << std::endl;
        std::cout << "  Growth D1(z)/D1(0) = " << growthFactor << std::endl;
    }
    
    // Build the inhomogeneous cosmology we will use.
    cosmo::PowerSpectrumPtr power;
    if(0 == bbandCoef) {
        // Use Eisenstein & Hu 1997 parameterization.
        boost::shared_ptr<cosmo::BaryonPerturbations> baryonsPtr(
            new cosmo::BaryonPerturbations(
            OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,baoOption));
        if(verbose) {
            // Print inhomogeneous cosmology info.
            std::cout << "z(eq) = " << baryonsPtr->getMatterRadiationEqualityRedshift() << std::endl;
            std::cout << "k(eq) = " << baryonsPtr->getMatterRadiationEqualityScale() << " /(Mpc/h)"
                << std::endl;
            std::cout << "sound horizon = " << baryonsPtr->getSoundHorizon() << " Mpc/h at z(drag) = "
                << baryonsPtr->getDragEpoch() << std::endl;
            std::cout << "Silk damping scale = " << baryonsPtr->getSilkDampingScale() << " /(Mpc/h)"
                << std::endl;
            double Tfc,Tfb,Tf;
            baryonsPtr->calculateTransferFunctions(kval,Tfc,Tfb,Tf);
            std::cout << "At k = " << kval << " /(Mpc/h):" << std::endl;
            std::cout << "  Tf(CDM,k) = " << Tfc << std::endl;
            std::cout << "  Tf(baryon,k) = " << Tfb << std::endl;
            std::cout << "  Tf(full,k) = " << Tf << std::endl;
        }

        // Create a sharable pointer to the matter transfer function (this will
        // keep baryonsPtr alive)
        cosmo::TransferFunctionPtr transferPtr(new cosmo::TransferFunction(boost::bind(
            &cosmo::BaryonPerturbations::getMatterTransfer,baryonsPtr,_1)));

        // Use COBE  n=1 normalization by default
        double deltaH = 1.94e-5*std::pow(OmegaMatter,-0.785-0.05*std::log(OmegaMatter));

        // Create a sharable pointer to a power spectrum for this transfer function
        // (this will keep transferPtr and therefore also baryonsPtr alive)
        boost::shared_ptr<cosmo::TransferFunctionPowerSpectrum> transferPowerPtr(
            new cosmo::TransferFunctionPowerSpectrum(transferPtr,spectralIndex,deltaH));

        // Use the requested sigma8 value if there is one.
        if(sigma8 > 0) {
            if(verbose) std::cout << "Renormalizing to sigma8 = " << sigma8 << std::endl;
            transferPowerPtr->setSigma(sigma8);
        }
        else if(verbose) {
            std::cout << "Using COBE n=1 deltaH = " << deltaH << std::endl;            
        }
        
        // Remember this power spectrum (this will keep all of the above alive)
        power = createFunctionPtr(transferPowerPtr);

        if(verbose) {
            std::cout << "Calculated sigma8(z=0) = " << cosmo::getRmsAmplitude(power,8)
                << std::endl;
        }

        // Rescale the power from z=0 to the desired redshift.
        transferPowerPtr->setDeltaH(growthFactor*transferPowerPtr->getDeltaH());

        if(0 < loadPowerFile.length()) {
            std::vector<std::vector<double> > columns(2);
            std::ifstream in(loadPowerFile.c_str());
            lk::readVectors(in,columns);
            in.close();
            if(verbose) {
                std::cout << "Read " << columns[0].size() << " rows from " << loadPowerFile
                    << std::endl;
            }
            lk::InterpolatorPtr iptr(new lk::Interpolator(columns[0],columns[1],"cspline"));
            power = cosmo::createFunctionPtr(iptr);
        }
    }
    else {
        if(0 < loadPowerFile.length()) {
            std::cerr << "Cannot combine --load-power with broadband options." << std::endl;
            return -1;
        }
        // Use a broadband power model.
        boost::shared_ptr<cosmo::BroadbandPower> bbPowerPtr(new cosmo::BroadbandPower(
            bbandCoef,bbandP,bbandKmin,bbandRmin,bbandR0,bbandVar));
        power = cosmo::createFunctionPtr(bbPowerPtr);
    }

    if(0 < savePowerFile.length()) {
        double pi(4*std::atan(1)),twopi2(2*pi*pi);
        boost::shared_ptr<cosmo::OneDimensionalPowerSpectrum> onedZero,onedHard,onedSoft;
        if(power1d) {
            onedZero.reset(new cosmo::OneDimensionalPowerSpectrum(power,0,kmin,kmax,nk));
            onedHard.reset(new cosmo::OneDimensionalPowerSpectrum(power,+r1d,kmin,kmax,nk));
            onedSoft.reset(new cosmo::OneDimensionalPowerSpectrum(power,-r1d,kmin,kmax,nk));
        }
        std::ofstream out(savePowerFile.c_str());
        double kratio(std::pow(kmax/kmin,1/(nk-1.)));
        for(int i = 0; i < nk; ++i) {
            double k(kmin*std::pow(kratio,i));
            if(k > kmax) k = kmax; // might happen with rounding
            out << k << ' ' << twopi2/(k*k*k)*(*power)(k);
            if(power1d) {
                out << ' ' << pi/k*(*onedZero)(k) << ' ' << pi/k*(*onedHard)(k)
                    << ' ' << pi/k*(*onedSoft)(k) << std::endl;
            }
        }
        out.close();
        if(verbose) {
            std::cout << "Wrote power to " << savePowerFile << std::endl;
        }
    }
    
    if(0 < saveCorrelationFile.length()) {
        cosmo::PowerSpectrumCorrelationFunction xi(power,rmin,rmax,multipole,nr);
        std::ofstream out(saveCorrelationFile.c_str());
        double r,dr;
        dr = rlog ? std::pow(rmax/rmin,1/(nr-1.)) : (rmax-rmin)/(nr-1.);
        for(int i = 0; i < nr; ++i) {
            r = rlog ? rmin*std::pow(dr,i) : rmin + dr*i;
            if(r > rmax) r = rmax; // might happen with rounding but xi(r) will complain
            out << r << ' ' << xi(r) << std::endl;
        }
        out.close();        
        if(verbose) {
            std::cout << "Wrote correlation function to " << saveCorrelationFile << std::endl;
        }
    }

    return 0;
}