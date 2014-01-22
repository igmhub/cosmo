// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A driver and test program for the MultipoleTransform class.
// Calculates the 3D or 2D transform of a tabulated multipole read from a file.

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
    po::options_description cli("Cosmology multipole transforms");
    std::string input,output;
    int ell,minSamplesPerCycle,minSamplesPerDecade;
    double min,max,veps,maxRelError;
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
        ("veps", po::value<double>(&veps)->default_value(1e-3),
            "desired transform accuracy")
        ("measure", "does initial measurements to optimize FFT plan")
        ("min-samples-per-cycle", po::value<int>(&minSamplesPerCycle)->default_value(2),
            "minimum number of samples per cycle to use for transform convolution")
        ("min-samples-per-decade", po::value<int>(&minSamplesPerDecade)->default_value(40),
            "minimum number of samples per decade to use for transform convolution")
        ("max-rel-error", po::value<double>(&maxRelError)->default_value(1e-3),
            "maximum allowed relative error for power-law extrapolation of input P(k)")
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
        measure(vm.count("measure"));

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

    cosmo::MultipoleTransform::Strategy strategy(measure ?
        cosmo::MultipoleTransform::MeasurePlan :
        cosmo::MultipoleTransform::EstimatePlan);

    try {
    	cosmo::MultipoleTransform mt(ttype,ell,min,max,veps,strategy,
            minSamplesPerCycle,minSamplesPerDecade);
        std::vector<double> const& ugrid = mt.getUGrid(), vgrid = mt.getVGrid();
        if(verbose) {
            std::cout << "Truncation fraction is " << mt.getTruncationFraction() << std::endl;
            std::cout << "Transform evaluated at " << mt.getNumPoints() << " points." << std::endl;
            std::cout <<  "Will evaluate at " << ugrid.size() << " points covering "
                << ugrid.front() << " to " << ugrid.back() << std::endl;
            std::cout <<  "Results estimated at " << vgrid.size() << " points covering "
                << vgrid.front() << " to " << vgrid.back() << std::endl;
        }
        std::vector<double> funcData(ugrid.size());
        for(int i = 0; i < funcData.size(); ++i) {
            funcData[i] = (*PkPtr)(ugrid[i]);
        }
        std::vector<double> results(vgrid.size());
        mt.transform(funcData,results);
        if(output.length() > 0) {
            std::ofstream out(output.c_str());
            for(int i = 0; i < results.size(); ++i) {
                out << vgrid[i] << ' ' << results[i] << std::endl;
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