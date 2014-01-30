// Created 22-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A driver and test program for the DistortedPowerCorrelation class.
// Calculates the 3D correlation function corresponding to a distorted power spectrum.

// Sample timing results using:
//
// ./cosmodpc -i ../../baofit/models/PlanckWPBestFitLCDM_matterpower.dat OPTIONS
//
//   TIME  OPTIONS
//  90.26s --repeat 1000
//  94.50s --repeat 1000 --optimize
//  64.13s --repeat 1000 --bypass
//  64.52s --repeat 1000 --bypass --optimize

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

class AutoCorrelationDistortion {
// A simple distortion model for autocorrelations, including linear redshift space effects
// (bias,beta), non-linear large-scale broadening (snlPar,snlPerp) and a continuum fitting
// broadband distortion model (k0,sigk).
public:
    AutoCorrelationDistortion(double bias, double beta,
        double snlPar, double snlPerp, double k0, double sigk) :
        _bias(bias), _beta(beta), _snlPar2(snlPar*snlPar), _snlPerp2(snlPerp*snlPerp),
        _k0(k0), _sigk(sigk)
    {
        _distScale = sigk > 0 ? 1/(1 + std::tanh(k0/sigk)) : 0;
    }
    double operator()(double k, double mu) const {
        double mu2(mu*mu);
        double linear = _bias*(1 + _beta*mu2);
        double snl2 = _snlPar2*mu2 + (1 - mu2)*_snlPerp2;
        double nonlinear = std::exp(-0.5*k*k*snl2);
        double kpar = std::fabs(k*mu);
        double distortion = 1 - _distScale*(1 - std::tanh((kpar-_k0)/_sigk));
        return distortion*nonlinear*linear*linear;
    }
private:
    double _bias,_beta,_snlPar2,_snlPerp2,_k0,_sigk,_distScale;
};

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology distorted power correlation function");
    std::string input,delta,output;
    int ellMax,nr,repeat,nk,nmu,minSamplesPerDecade;
    double rmin,rmax,relerr,abserr,abspow,maxRelError,kmin,kmax,margin,vepsMin,vepsMax;
    double bias,beta,snlPar,snlPerp,k0,sigk;
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
        ("bias", po::value<double>(&bias)->default_value(1.0),
            "linear tracer bias")
        ("beta", po::value<double>(&beta)->default_value(1.4),
            "redshift-space distortion parameter")
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
        optimize(vm.count("optimize")), bypass(vm.count("bypass"));

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

        boost::shared_ptr<AutoCorrelationDistortion> rsd(
            new AutoCorrelationDistortion(bias,beta,snlPar,snlPerp,k0,sigk));
        cosmo::RMuFunctionCPtr distPtr(new cosmo::RMuFunction(boost::bind(
            &AutoCorrelationDistortion::operator(),rsd,_1,_2)));

    	cosmo::DistortedPowerCorrelation dpc(PkPtr,distPtr,rmin,rmax,nr,ellMax,
            symmetric,relerr,abserr,abspow);
        // initialize
        dpc.initialize(nmu,minSamplesPerDecade,margin,vepsMax,vepsMin,optimize);
        if(verbose) dpc.printToStream(std::cout);
        // transform (with repeats, if requested)
        bool ok;
        for(int i = 0; i < repeat; ++i) {
            ok = dpc.transform(bypass);
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
                    kout << ' ' << boost::lexical_cast<std::string>(dpc.getPowerMultipole(k,ell));
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