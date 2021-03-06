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
#include <vector>
#include <list>

namespace po = boost::program_options;
namespace lk = likely;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology calculator");
    double OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,spectralIndex,sigma8,
        zval,kval,kmin,kmax,r1d,rmin,rmax,baoAmplitude,baoSigma,baoScale;
    double bbandP,bbandCoef,bbandKmin,bbandRmin,bbandR0,bbandVar,epsAbs,epsRel;
    int nk,nr;
    std::string loadPowerFile,savePowerFile,saveCorrelationFile;
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0.272),
            "Present-day value of OmegaMatter.")
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
        ("bao-smooth", "Smooths out any BAO oscillation by sampling only at nodes.")
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
        ("eps-abs", po::value<double>(&epsAbs)->default_value(1e-8),
            "Absolute precision goal for 1D integration of correlation functions.")
        ("eps-rel", po::value<double>(&epsRel)->default_value(1e-6),
            "Relative precision goal for 1D integration of correlation functions (integrand1 only).")
        ("rlog", "Use log spaced r-values for saved correlation function (default is linear).")
        ("quad", "Calculates the quadrupole (l=2) correlation function (default is monopole).")
        ("hexa", "Calculates the hexedacapole (l=4) correlation function (default is monopole).")
        ("no-wiggles", "Calculates the power spectrum without baryon acoustic oscillations.")
        ("no-osc", "Calculates the power spectrum with a non-oscillation baryon transfer function.")
        ("periodic-osc", "Calculates the power spectrum with periodic acoustic oscillations.")
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
        noOsc(vm.count("no-osc")),periodicOsc(vm.count("periodic-osc")), baoSmooth(vm.count("bao-smooth"));

    // Process the multipole flags.
    if(quad && hexa) {
        std::cerr << "Cannot request both quad (l=2) and hexa (l=4) for correlation function output."
            << std::endl;
        return -1;
    }
    cosmo::Multipole multipole(cosmo::Monopole);
    if(quad) multipole = cosmo::Quadrupole;
    if(hexa) multipole = cosmo::Hexadecapole;

    // Process the wiggle flags.
    if(vm.count("no-wiggles")+vm.count("periodic-wiggles")+vm.count("bao-fit") > 1) {
        std::cerr << "Specify at most one of no-wiggles, periodic-wiggles, bao-fit options."
            << std::endl;
        return -1;
    }
    cosmo::BaryonPerturbations::BaoOption baoOption(cosmo::BaryonPerturbations::ShiftedOscillation);
    if(noOsc) baoOption = cosmo::BaryonPerturbations::NoOscillation;
    if(periodicOsc) baoOption = cosmo::BaryonPerturbations::PeriodicOscillation;

    try {
        // Build a homogeneous cosmology using the parameters specified (and OmegaK = 0).
        cosmo::LambdaCdmRadiationUniverse cosmology(OmegaMatter,0,hubbleConstant,cmbTemp);
        double growthFactor = cosmology.getGrowthFunction(zval)/cosmology.getGrowthFunction(0);
        if(verbose) {
            // Print homogeneous cosmology info.
            std::cout << "OmegaRadiation = " << cosmology.getOmegaRadiation() << std::endl;
            std::cout << "OmegaLambda = " << cosmology.getOmegaLambda() << std::endl;
            std::cout << "OmegaK = " << cosmology.getCurvature() << std::endl;    
            std::cout << "At z = " << zval << ':' << std::endl;
            std::cout << "  H(z)/H0 = " << cosmology.getHubbleFunction(zval) << std::endl;
            std::cout << "  Radial D(z) = " << cosmology.getLineOfSightComovingDistance(zval)
                << " Mpc/h" << std::endl;
            std::cout << "  Transverse DA(z) = " << cosmology.getTransverseComovingScale(zval)
                << " Mpc/h/rad" << std::endl;
            double tL(cosmology.getLookbackTime(zval));
            double conv(1e9*86400*365.25);
            std::cout << "  t(lookback,z) = " << tL << " secs/h = " << tL/conv*hubbleConstant
                << " Gyr/h" << std::endl;
            std::cout << "  Growth D1(z)/D1(0) = " << growthFactor << std::endl;
        }
    
        // Build the inhomogeneous cosmology we will use.
        cosmo::PowerSpectrumPtr power;
        if(0 == bbandCoef) {
            // Create an Eisenstein & Hu 1997 parameterization.
            boost::shared_ptr<cosmo::BaryonPerturbations> baryonsPtr(
                new cosmo::BaryonPerturbations(
                OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,baoOption));
            if(verbose) {
                // Print inhomogeneous cosmology info.
                std::cout << "z(eq) = " << baryonsPtr->getMatterRadiationEqualityRedshift() << std::endl;
                std::cout << "k(eq) = " << baryonsPtr->getMatterRadiationEqualityScale() << " /(Mpc/h)"
                    << std::endl;
                std::cout << "sound horizon = " << baryonsPtr->getSoundHorizon() << " Mpc/h at z(drag) = "
                    << baryonsPtr->getDragEpoch() << " (fit = " << baryonsPtr->getSoundHorizonFit()
                    << " Mpc/h)" << std::endl;
                std::cout << "Silk damping scale = " << baryonsPtr->getSilkDampingScale() << " /(Mpc/h)"
                    << std::endl;
                double Tfc,Tfb,Tf,Tnw;
                baryonsPtr->calculateTransferFunctions(kval,Tfb,Tfc,Tf,Tnw,baoOption);
                std::cout << "At k = " << kval << " /(Mpc/h):" << std::endl;
                std::cout << "  Tf(CDM,k) = " << Tfc << std::endl;
                std::cout << "  Tf(baryon,k) = " << Tfb << std::endl;
                std::cout << "  Tf(full,k) = " << Tf << std::endl;
                std::cout << "  Tf(nw,k) = " << Tnw << std::endl;
            }

            if(0 < loadPowerFile.length()) {
                // Load a tabulated power spectrum for interpolation.
                std::vector<std::vector<double> > columns(2);
                std::ifstream in(loadPowerFile.c_str());
                lk::readVectors(in,columns);
                in.close();
                if(verbose) {
                    std::cout << "Read " << columns[0].size() << " rows from " << loadPowerFile
                        << std::endl;
                }
                double pi(4*std::atan(1)),twopi2(2*pi*pi);
                // rescale to k^3/(2pi^2) P(k)
                for(int row = 0; row < columns[0].size(); ++row) {
                    double k(columns[0][row]);
                    columns[1][row] *= k*k*k/twopi2;
                }
                // Create an interpolator of this data.
                lk::InterpolatorPtr iptr(new lk::Interpolator(columns[0],columns[1],"cspline"));
                // Use the resulting interpolation function for future power calculations.
                power = likely::createFunctionPtr(iptr);
            }
            else {
                // Create a power spectrum from the EH97 parameterization...

                // Create a sharable pointer to the matter transfer function (this will
                // keep baryonsPtr alive)
                cosmo::TransferFunctionPtr transferPtr(new cosmo::TransferFunction(boost::bind(
                    noWiggles ?
                        &cosmo::BaryonPerturbations::getNoWigglesTransfer :
                        &cosmo::BaryonPerturbations::getMatterTransfer,
                    baryonsPtr,_1)));

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
                power = likely::createFunctionPtr(transferPowerPtr);

                if(verbose) {
                    std::cout << "Calculated sigma8(z=0) = " << cosmo::getRmsAmplitude(power,8)
                        << std::endl;
                }

                // Rescale the power from z=0 to the desired redshift.
                transferPowerPtr->setDeltaH(growthFactor*transferPowerPtr->getDeltaH());
            }

            if(verbose) {
                std::cout << "Calculated sigma8(z=" << zval << ") = " << cosmo::getRmsAmplitude(power,8)
                    << std::endl;
            }
        
            // Resample at BAO nodes if requested
            if(baoSmooth) {
                std::list<double> klist;
                double pi(4*std::atan(1));
                double ks = pi/baryonsPtr->getSoundHorizon();
                double k1 = ks*baryonsPtr->getNode(1);
                if(kmax < k1) {
                    std::cerr << "Cannot smooth with kmax < k1.";
                    return -2;
                }
                // Sample at the first 20 nodes.
                int numNodes(20);
                for(int n = 1; n < numNodes; ++n) {
                    double k = ks*baryonsPtr->getNode(n);
                    klist.push_back(k);
                }
                // Sample up to kmax doubling the node number between samples.
                double k;
                int n(numNodes);
                double smoothStep(1.25);
                do {
                    k = ks*baryonsPtr->getNode(n);
                    klist.push_back(k);
                    n = (int)std::ceil(n*smoothStep);
                } while(k < kmax);
                // Sample from k1 down to kmin, halving k each time.
                k = k1;
                do {
                    k /= smoothStep;
                    klist.push_front(k);
                } while(k > kmin);
                // Calculate the power at each k value.
                std::list<double>::iterator nextk(klist.begin());
                std::vector<double> kvec(klist.size()), pvec(klist.size());
                for(int i = 0; i < kvec.size(); ++i) {
                    k = kvec[i] = *nextk++;
                    pvec[i] = (*power)(kvec[i]);
                }
                // Create an interpolator.
                lk::InterpolatorPtr iptr(new lk::Interpolator(kvec,pvec,"cspline"));
                // Use the resulting interpolation function for future power calculations.
                power = likely::createFunctionPtr(iptr);
            }
        }
        else {
            // Create a broadband power model.
            if(0 < loadPowerFile.length()) {
                std::cerr << "Cannot combine --load-power with broadband options." << std::endl;
                return -1;
            }
            boost::shared_ptr<cosmo::BroadbandPower> bbPowerPtr(new cosmo::BroadbandPower(
                bbandCoef,bbandP,bbandKmin,bbandRmin,bbandR0,bbandVar));
            power = likely::createFunctionPtr(bbPowerPtr);
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
                        << ' ' << pi/k*(*onedSoft)(k);
                }
                out << std::endl;
            }
            out.close();
            if(verbose) {
                std::cout << "Wrote power to " << savePowerFile << std::endl;
            }
        }
    
        if(0 < saveCorrelationFile.length()) {
            cosmo::PowerSpectrumCorrelationFunction xi(power,rmin,rmax,multipole,nr,epsAbs,epsRel);
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
    }
    catch(std::runtime_error const &e) {
        std::cerr << "ERROR: exiting with an exception:\n  " << e.what() << std::endl;
        return -4;
    }

    return 0;
}
