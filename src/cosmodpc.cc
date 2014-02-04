// Created 22-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A driver and test program for the DistortedPowerCorrelation class.
// Calculates the 3D correlation function corresponding to a distorted power spectrum.

// Sample timing results on OS 10.8 laptop using:
//
// ./cosmodpc -i ../../baofit/models/PlanckWPBestFitLCDM_matterpower.dat OPTIONS
//
//   TIME  OPTIONS
//  127.3s --repeat 10000 --direct-power-multipoles
//  126.5s --repeat 10000 --direct-power-multipoles --optimize
//   83.2s --repeat 10000 --direct-power-multipoles --bypass
//   53.2s --repeat 10000
//   53.4s --repeat 10000 --optimize
//   53.2s --repeat 10000 --bypass

#include "cosmo/cosmo.h"
#include "likely/likely.h"
#include "likely/function_impl.h"

#include "boost/program_options.hpp"
#include "boost/bind.hpp"
#include "boost/lexical_cast.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

namespace po = boost::program_options;
namespace lk = likely;

class LyaDistortion {
// A simple distortion model for autocorrelations, including linear redshift space effects
// (bias,beta), non-linear large-scale broadening (snlPar,snlPerp) and a continuum fitting
// broadband distortion model (k0,sigk).
public:
    LyaDistortion(double bias, double biasbeta,
        double biasGamma, double biasSourceAbsorber, double biasAbsorberResponse,
        double meanFreePath, double snlPar, double snlPerp, double k0, double sigk) :
        _bias(bias), _biasbeta(biasbeta),
        _biasGamma(biasGamma), _biasSourceAbsorber(biasSourceAbsorber),
        _biasAbsorberResponse(biasAbsorberResponse), _meanFreePath(meanFreePath),
        _snlPar2(snlPar*snlPar), _snlPerp2(snlPerp*snlPerp),
        _k0(k0), _sigk(sigk)
    {
        _radStrength = biasGamma*biasSourceAbsorber;
        _distScale = sigk > 0 ? 1/(1 + std::tanh(k0/sigk)) : 0;
    }
    double operator()(double k, double mu) const {
        // Calculate the k-dependent effective bias
        double bias(_bias);
        if(_radStrength != 0) {
            double s(_meanFreePath*k);
            double Ws = std::atan(s)/s;
            bias += _radStrength*Ws/(1+_biasAbsorberResponse*Ws);
        }
        // Calculate the k-dependent effective beta
        double beta(_biasbeta/bias);
        // Calculate the overall large-scale Lya tracer bias
        double mu2(mu*mu);
        double linear = bias*(1 + beta*mu2);
        // Calculate non-linear broadening
        double snl2 = _snlPar2*mu2 + (1 - mu2)*_snlPerp2;
        double nonlinear = std::exp(-0.5*k*k*snl2);
        // Calculate continuum fitting distortion
        double kpar = std::fabs(k*mu);
        double distortion = 1 - _distScale*(1 - std::tanh((kpar-_k0)/_sigk));
        return distortion*nonlinear*linear*linear;
    }
private:
    double _bias,_biasbeta,_biasGamma,_biasSourceAbsorber,_biasAbsorberResponse,_meanFreePath,
        _snlPar2,_snlPerp2,_k0,_sigk,_distScale,_radStrength;
};

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology distorted power correlation function");
    std::string input,delta,output;
    int ellMax,nr,repeat,nk,nmu,minSamplesPerDecade;
    double rmin,rmax,relerr,abserr,abspow,maxRelError,kmin,kmax,margin,vepsMin,vepsMax;
    double bias,biasbeta,biasGamma,biasSourceAbsorber,biasAbsorberResponse,meanFreePath,
        snlPar,snlPerp,k0,sigk;
    cli.add_options()
        ("help,h", "prints this info and exits.")
        ("verbose", "prints additional information.")
        ("input,i", po::value<std::string>(&input)->default_value(""),
            "filename to read k,P(k) values from")
        ("delta", po::value<std::string>(&delta)->default_value(""),
            "optional filename of k,P(k) values to subtract from input")
        ("output,o", po::value<std::string>(&output)->default_value(""),
            "base name for saving results")
        ("rmin", po::value<double>(&rmin)->default_value(10.),
            "minimum value of comoving separation to use")
        ("rmax", po::value<double>(&rmax)->default_value(200.),
            "maximum value of comoving separation to use")
        ("nr", po::value<int>(&nr)->default_value(191),
            "number of points spanning [rmin,rmax] to use")
        ("ell-max", po::value<int>(&ellMax)->default_value(4),
            "maximum multipole to use for transforms")
        ("asymmetric", "distortion is asymmetric in mu (uses odd ell values)")
        ("bias", po::value<double>(&bias)->default_value(-0.17),
            "linear tracer bias")
        ("biasbeta", po::value<double>(&biasbeta)->default_value(-0.17),
            "product of bias and linear redshift-space distortion parameter beta")
        ("bias-gamma", po::value<double>(&biasGamma)->default_value(0),
            "Lyman-alpha photoionization bias")
        ("bias-source-absorber", po::value<double>(&biasSourceAbsorber)->default_value(1.0),
            "Lyman-alpha source - absorber bias difference")
        ("bias-absorber-response", po::value<double>(&biasAbsorberResponse)->default_value(-2./3.),
            "Lyman-alpha absorber response bias")
        ("mean-free-path", po::value<double>(&meanFreePath)->default_value(300.),
            "effective mean free path of an ionizing photon in Mpc/h")
        ("snl-par", po::value<double>(&snlPar)->default_value(0.),
            "parallel component of non-linear broadening in Mpc/h")
        ("snl-perp", po::value<double>(&snlPerp)->default_value(0.),
            "perpendicular component of non-linear broadening in Mpc/h")
        ("k0", po::value<double>(&k0)->default_value(0.02),
            "cutoff scale for broadband distortion in h/Mpc (ignored when sigk = 0)")
        ("sigk", po::value<double>(&sigk)->default_value(0.),
            "smoothing scale for broadband distortion in h/Mpc (or zero for no distortion)")
        ("relerr", po::value<double>(&relerr)->default_value(1e-3),
            "relative error termination goal")
        ("abserr", po::value<double>(&abserr)->default_value(1e-5),
            "absolute error termination goal")
        ("abspow", po::value<double>(&abspow)->default_value(0.),
            "absolute error weighting power")
        ("margin", po::value<double>(&margin)->default_value(2.),
            "termination criteria margin to use for initialization")
        ("veps-max", po::value<double>(&vepsMax)->default_value(0.01),
            "maximum value of veps value to use")
        ("veps-min", po::value<double>(&vepsMin)->default_value(1e-6),
            "minimum value of veps value to use")
        ("min-samples-per-decade", po::value<int>(&minSamplesPerDecade)->default_value(40),
            "minimum number of samples per decade to use for transform convolution")
        ("max-rel-error", po::value<double>(&maxRelError)->default_value(1e-3),
            "maximum allowed relative error for power-law extrapolation of input P(k)")
        ("direct-power-multipoles",
            "use direct calculation of P(k) multipoles instead of interpolation")
        ("optimize", "optimizes transform FFTs")
        ("bypass", "bypasses the termination test for transforms")
        ("repeat", po::value<int>(&repeat)->default_value(1),
            "number of times to repeat identical transform")
        ("kmin", po::value<double>(&kmin)->default_value(0.005),
            "minimum value of comoving separation to use")
        ("kmax", po::value<double>(&kmax)->default_value(0.35),
            "maximum value of comoving separation to use")
        ("nk", po::value<int>(&nk)->default_value(100),
            "number of log-spaced k values for saving results")
        ("nmu", po::value<int>(&nmu)->default_value(10),
            "number of equally spaced mu_k and mu_r values for saving results")
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
    bool verbose(vm.count("verbose")), symmetric(0==vm.count("asymmetric")),
        optimize(vm.count("optimize")), bypass(vm.count("bypass")),
        directPowerMultipoles(vm.count("direct-power-multipoles"));

    if(!symmetric) {
        std::cerr << "Odd multipoles not implemented yet." << std::endl;
        return 1;
    }

    if(input.length() == 0) {
        std::cerr << "Missing input filename." << std::endl;
        return 1;
    }

    int dell(symmetric ? 2:1);

    try {
        cosmo::TabulatedPowerCPtr power =
            cosmo::createTabulatedPower(input,true,true,maxRelError,verbose);
        if(delta.length() > 0) {
            cosmo::TabulatedPowerCPtr power2 =
                cosmo::createTabulatedPower(delta,true,true,maxRelError,verbose);
            power = power->createDelta(power2,verbose);
        }
        lk::GenericFunctionPtr PkPtr =
            lk::createFunctionPtr<const cosmo::TabulatedPower>(power);

        boost::shared_ptr<LyaDistortion> rsd(new LyaDistortion(
            bias,biasbeta,biasGamma,biasSourceAbsorber,biasAbsorberResponse,meanFreePath,
            snlPar,snlPerp,k0,sigk));
        cosmo::RMuFunctionCPtr distPtr(new cosmo::RMuFunction(boost::bind(
            &LyaDistortion::operator(),rsd,_1,_2)));

        // Use the limits of the input tabulated power for tabulating the
        // power multipoles (the kmin,kmax cmd-line args are for output only)
        double klo = power->getKMin(), khi = power->getKMax();
        int nkint = std::ceil(std::log10(khi/klo)*minSamplesPerDecade);
    	cosmo::DistortedPowerCorrelation dpc(PkPtr,distPtr,
            klo,khi,nkint,rmin,rmax,nr,ellMax,
            symmetric,relerr,abserr,abspow);
        // initialize
        dpc.initialize(nmu,margin,vepsMax,vepsMin,optimize);
        if(verbose) dpc.printToStream(std::cout);
        // transform (with repeats, if requested)
        bool ok;
        for(int i = 0; i < repeat; ++i) {
            ok = dpc.transform(!directPowerMultipoles,bypass);
        }
        if(!ok) {
            std::cerr << "Transform fails termination test." << std::endl;
        }
        // save transform results
        if(output.length() > 0) {
            int dell = symmetric ? 2 : 1;
            double dmu = symmetric ? 1./(nmu-1.) : 2./(nmu-1.);
            // Write out values tabulated for log-spaced k
            double dk = std::pow(kmax/kmin,1./(nk-1.));
            std::string kfile = output + ".k.dat";
            std::ofstream kout(kfile.c_str());
            for(int i = 0; i < nk; ++i) {
                double k = kmin*std::pow(dk,i);
                kout << boost::lexical_cast<std::string>(k);
                for(int j = 0; j < nmu; ++j) {
                    double mu = 1 - j*dmu;
                    kout << ' ' << boost::lexical_cast<std::string>(dpc.getPower(k,mu));
                }
                for(int ell = 0; ell <= ellMax; ell += dell) {
                    double pk = (directPowerMultipoles) ?
                        dpc.getPowerMultipole(k,ell) : dpc.getSavedPowerMultipole(k,ell);
                    kout << ' ' << boost::lexical_cast<std::string>(pk);
                }
                kout << std::endl;
            }
            kout.close();
            // Write out values tabulated for linear-spaced r
            double dr = (rmax-rmin)/(nr-1.);
            std::string rfile = output + ".r.dat";
            std::ofstream rout(rfile.c_str());
            for(int i = 0; i < nr; ++i) {
                double r = rmin + i*dr;
                rout << boost::lexical_cast<std::string>(r);
                for(int j = 0; j < nmu; ++j) {
                    double mu = 1 - j*dmu;
                    rout << ' ' << boost::lexical_cast<std::string>(dpc.getCorrelation(r,mu));
                }
                for(int ell = 0; ell <= ellMax; ell += dell) {
                    rout << ' ' << boost::lexical_cast<std::string>(dpc.getCorrelationMultipole(r,ell));
                }
                rout << std::endl;                
            }
            rout.close();
        }
    }
    catch(std::runtime_error const &e) {
        std::cerr << "ERROR: exiting with an exception:\n  " << e.what() << std::endl;
        return -1;
    }

    return 0;
}