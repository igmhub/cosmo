// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A driver and test program for the AdaptiveMultipoleTransform class.
// Calculates the 3D or 2D transform of a tabulated multipole read from a file.

// Sample timing results using:
//
// ./cosmoatrans -i ../../baofit/models/PlanckWPBestFitLCDM_matterpower.dat OPTIONS
//
//   TIME  OPTIONS
//  15.26s --repeat 100000
//  11.16s --repeat 100000 --optimize
//   9.82s --repeat 100000 --bypass
//   6.69s --repeat 100000 --optimize --bypass

#include "cosmo/cosmo.h"
#include "likely/likely.h"
#include "likely/function_impl.h"

#include "boost/program_options.hpp"
#include "boost/bind.hpp"

#include <iostream>
#include <fstream>

namespace po = boost::program_options;
namespace lk = likely;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology adaptive multipole transforms");
    std::string input,output;
    int ell,npoints,minSamplesPerDecade,repeat;
    double min,max,relerr,abserr,abspow,margin,vepsMax,vepsMin,maxRelError;
    cli.add_options()
        ("help,h", "prints this info and exits.")
        ("verbose", "prints additional information.")
        ("input,i", po::value<std::string>(&input)->default_value(""),
            "name of filename to read k,P_ell(k) values from")
        ("output,o", po::value<std::string>(&output)->default_value(""),
            "name of filename to save r,xi_ell(r) values to")
        ("hankel", "performs a 2D Hankel transform (default is 3D spherical Bessel")
        ("ell", po::value<int>(&ell)->default_value(0),
            "multipole number of transform to calculate")
        ("min", po::value<double>(&min)->default_value(10.),
            "minimum value of transformed coordinate")
        ("max", po::value<double>(&max)->default_value(200.),
            "maximum value of transformed coordinate")
        ("npoints", po::value<int>(&npoints)->default_value(50),
            "number of points spanning [min,max] to use")
        ("relerr", po::value<double>(&relerr)->default_value(1e-2),
            "relative error termination goal")
        ("abserr", po::value<double>(&abserr)->default_value(1e-3),
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
    bool verbose(vm.count("verbose")),hankel(vm.count("hankel")),
        optimize(vm.count("optimize")), bypass(vm.count("bypass"));

    if(input.length() == 0) {
        std::cerr << "Missing input filename." << std::endl;
        return 1;
    }
    cosmo::TabulatedPowerCPtr power =
        cosmo::createTabulatedPower(input,true,true,maxRelError,verbose);
    lk::GenericFunctionPtr PkPtr =
        lk::createFunctionPtr<const cosmo::TabulatedPower>(power);

    cosmo::MultipoleTransform::Type ttype(hankel ?
        cosmo::MultipoleTransform::Hankel :
        cosmo::MultipoleTransform::SphericalBessel);

    std::vector<double> points;
    double dv = (max-min)/(npoints-1.);
    for(int i = 0; i < npoints; ++i) {
        points.push_back(min + i*dv);
    }

    try {
    	cosmo::AdaptiveMultipoleTransform mt(ttype,ell,points,relerr,abserr,abspow);
        std::vector<double> result(npoints);
        double veps = mt.initialize(PkPtr,result,minSamplesPerDecade,margin,
            vepsMax,vepsMin,optimize);
        if(verbose) {
            std::cout << "Using veps = " << veps << std::endl;
        }
        bool ok;
        for(int i = 0; i < repeat; ++i) {
            ok = mt.transform(PkPtr,result,bypass);
        }
        if(!ok) {
            std::cerr << "Transform fails termination test." << std::endl;
        }
        if(output.length() > 0) {
            std::ofstream out(output.c_str());
            for(int i = 0; i < npoints; ++i) {
                out << points[i] << ' ' << result[i] << std::endl;
            }
            out.close();
        }
    }
    catch(std::runtime_error const &e) {
        std::cerr << "ERROR: exiting with an exception:\n  " << e.what() << std::endl;
        return -1;
    }

    return 0;
}