// Created 09-Jul-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Generates a Gaussian random field.

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <iostream>
#include <fstream>
#include <string>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    double spacing;
    int nx,ny,nz,seed;
    std::string loadPowerFile;
    po::options_description cli("Gaussian random field generator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("spacing", po::value<double>(&spacing)->default_value(1),
            "Grid spacing in Mpc/h.")
        ("nx", po::value<int>(&nx)->default_value(64),
            "Grid size along x-axis.")
        ("ny", po::value<int>(&ny)->default_value(0),
            "Grid size along y-axis (or zero for ny=nx).")
        ("nz", po::value<int>(&nz)->default_value(0),
            "Grid size along z-axis (or zero for nz=ny).")
        ("load-power", po::value<std::string>(&loadPowerFile)->default_value(""),
            "Reads k,P(k) values (in h/Mpc units) to interpolate from the specified filename.")
        ("seed", po::value<int>(&seed)->default_value(123),
            "Random seed to use.")
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

    // Fill in any missing grid dimensions.
    if(0 == ny) ny = nx;
    if(0 == nz) nz = ny;
    
    double pi(4*std::atan(1));
    
    // Load a tabulated power spectrum for interpolation.
    cosmo::PowerSpectrumPtr power;
    if(0 < loadPowerFile.length()) {
        std::vector<std::vector<double> > columns(2);
        std::ifstream in(loadPowerFile.c_str());
        likely::readVectors(in,columns);
        in.close();
        if(verbose) {
            std::cout << "Read " << columns[0].size() << " rows from " << loadPowerFile
                << std::endl;
        }
        double twopi2(2*pi*pi);
        // rescale to k^3/(2pi^2) P(k)
        for(int row = 0; row < columns[0].size(); ++row) {
            double k(columns[0][row]);
            columns[1][row] *= k*k*k/twopi2;
        }
        // Create an interpolator of this data.
        likely::InterpolatorPtr iptr(new likely::Interpolator(columns[0],columns[1],"cspline"));
        // Use the resulting interpolation function for future power calculations.
        power = likely::createFunctionPtr(iptr);
    }
    else {
        std::cerr << "Missing required load-power filename." << std::endl;
        return -2;
    }
    
    // Initialize the random number source.
    likely::Random::instance()->setSeed(seed);
    
    // Create the generator.
    cosmo::FftGaussianRandomFieldGenerator generator(power, spacing, nx, ny, nz);
    if(verbose) {
        std::cout << "Memory size = "
            << boost::format("%.1f Mb") % (generator.getMemorySize()/1048576.) << std::endl;
    }
    generator.generate();
    
    // Calculate the statistics of the generated delta field.
    likely::WeightedAccumulator accumulator;
    for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
            for(int iz = 0; iz < nz; ++iz) {
                accumulator.accumulate(generator.getField(ix,iy,iz));
            }
        }
    }
    
    // Compare with var(k1,k2) = Integral[k^2/(2pi) P(k),{k,k1,k2}]. This will not match exactly
    // because we are comparing a sphere in k-space with a cuboid in r-space.
    double kmax = pi/spacing, kmin = pi/(spacing*std::pow(nx*ny*nz,1./3.));
    int nsteps(1000);
    double kratio = std::pow(kmax/kmin,1./(nsteps-1.));
    double kVariance(0);
    for(int step = 0; step < nsteps; ++step) {
        double k = kmin*std::pow(kratio,step);
        kVariance += (*power)(k)/k;
    }
    kVariance *= std::log(kratio);
    
    std::cout << "Delta field mean = " << accumulator.mean() << ", variance = "
        << accumulator.variance() << ", P(k) variance estimate is " << kVariance << std::endl;
    
    return 0;
}
