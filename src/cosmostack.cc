// Created 8-Oct-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Stacks many Gaussian random fields on the field maximum (or minimum).

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <iostream>
#include <fstream>
#include <string>

namespace po = boost::program_options;
namespace lk = likely;

// Return the distance between x0 and x1 in a periodic dimension of length nx
int distance(int x0, int x1, int nx){
    int dx(x1 - x0);
    return (dx > 0 ? (dx > nx/2 ? nx - dx : dx) : (dx <= -nx/2 ? nx + dx : dx) );
}

int main(int argc, char **argv) {
    // Configure command-line option processing
    double spacing, xlos, ylos, zlos;
    long npairs;
    int nx,ny,nz,seed,nfields;
    std::string loadPowerFile, prefix;
    po::options_description cli("Stacks many Gaussian random fields on the field maximum (or minimum).");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("spacing", po::value<double>(&spacing)->default_value(4),
            "Grid spacing in Mpc/h.")
        ("nx", po::value<int>(&nx)->default_value(76),
            "Grid size along x-axis.")
        ("ny", po::value<int>(&ny)->default_value(0),
            "Grid size along y-axis (or zero for ny=nx).")
        ("nz", po::value<int>(&nz)->default_value(0),
            "Grid size along z-axis (or zero for nz=ny).")
        ("load-power", po::value<std::string>(&loadPowerFile)->default_value(""),
            "Reads k,P(k) values (in h/Mpc units) to interpolate from the specified filename.")
        ("seed", po::value<int>(&seed)->default_value(511),
            "Random seed to use for GRF.")
        ("prefix", po::value<std::string>(&prefix)->default_value("stack"),
            "Prefix for output file names.")
        ("nfields", po::value<int>(&nfields)->default_value(1000),
            "Number of fields to stack.")
        ("fiducial",
            "Stack relative to box center instead of local maximum.")
        ("minimum",
            "Stack relative to the local minimum instead of local maximum.")
        ("snapshot",
            "Save snapshot of 1d projection every 10%%.")
        ("xlos", po::value<double>(&xlos)->default_value(1),
            "Line-of-sight direction x-component.")
        ("ylos", po::value<double>(&ylos)->default_value(0),
            "Line-of-sight direction y-component.")
        ("zlos", po::value<double>(&zlos)->default_value(0),
            "Line-of-sight direction z-component.")
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
    bool verbose(vm.count("verbose")), fiducial(vm.count("fiducial")), snapshot(vm.count("snapshot")),
        minimum(vm.count("minimum"));

    // Fill in any missing grid dimensions.
    if(0 == ny) ny = nx;
    if(0 == nz) nz = ny;
    
    double pi(4*std::atan(1));

    // Load a tabulated power spectrum for interpolation.
    cosmo::PowerSpectrumPtr power;
    if(0 < loadPowerFile.length()) {
        std::vector<std::vector<double> > columns(2);
        std::ifstream in(loadPowerFile.c_str());
        lk::readVectors(in,columns);
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
        lk::InterpolatorPtr iptr(new lk::Interpolator(columns[0],columns[1],"cspline"));
        // Use the resulting interpolation function for future power calculations.
        power = lk::createFunctionPtr(iptr);
    }
    else {
        std::cerr << "Missing required load-power filename." << std::endl;
        return -2;
    }

    // Initialize the random number source.
    lk::Random::instance()->setSeed(seed);

    // Create the generator.
    cosmo::FftGaussianRandomFieldGenerator generator(power, spacing, nx, ny, nz);
    if(verbose) {
        std::cout << "Memory size = "
            << boost::format("%.1f Mb") % (generator.getMemorySize()/1048576.) << std::endl;
    }

    // Initialize histogram
    double binspacing(spacing);
    int nbins(std::floor(150./binspacing + .5));
    double rmax(nbins*binspacing);
    boost::scoped_array<lk::WeightedAccumulator> xi(new lk::WeightedAccumulator[nbins]);
    boost::scoped_array<lk::WeightedAccumulator> xi2d(new lk::WeightedAccumulator[nbins*nbins]);

    int dx,dy,dz;
    double value, extremeValue, dr, r, rperp, rparl;
    std::vector<int> extremeIndex(3,0);
    lk::WeightedAccumulator extremeValues;
    double normlos(std::sqrt(xlos*xlos + ylos*ylos + zlos*zlos));
    double xparl(xlos/normlos), yparl(ylos/normlos), zparl(zlos/normlos);
    for(int ifield = 0; ifield < nfields; ++ifield){
        // Generate Gaussian random field
        generator.generate();
        extremeValue = generator.getField(0,0,0);
        if(!fiducial) {
            for(int iz = 0; iz < nz; ++iz){
                for(int iy = 0; iy < ny; ++iy){
                    for(int ix = 0; ix < nx; ix++){
                        value = generator.getField(ix,iy,iz);
                        if((minimum ? value < extremeValue : value > extremeValue)){
                            extremeValue = value;
                            extremeIndex[0] = ix;
                            extremeIndex[1] = iy;
                            extremeIndex[2] = iz;
                        }
                    }
                }
            }
        }
        else {
            extremeIndex[0] = nx/2; extremeIndex[1] = ny/2; extremeIndex[2] = nz/2;
        }
        // Accumulate extreme value
        extremeValues.accumulate(generator.getField(extremeIndex[0],extremeIndex[1],extremeIndex[2]));
        // Fill 1-d, 2-d histograms
        for(int iz = 0; iz < nz; ++iz){
            for(int iy = 0; iy < ny; ++iy){
                for(int ix = 0; ix < nx; ix++){
                    // Calculate distance from extreme value point on grid (wrap-around)
                    dx = distance(ix,extremeIndex[0],nx);
                    dy = distance(iy,extremeIndex[1],ny);
                    dz = distance(iz,extremeIndex[2],nz);
                    dr = std::sqrt(dx*dx + dy*dy + dz*dz);
                    // Apply grid spacing
                    r = spacing*dr;
                    rparl = spacing*std::abs(dx*xparl + dy*yparl + dz*zparl);
                    rperp = std::sqrt(r*r-rparl*rparl);
                    // Look up field value
                    value = generator.getField(ix,iy,iz);
                    // Accumulate value in appropriate bin
                    if(r < rmax) {
                        xi[std::floor(r/binspacing)].accumulate(value);
                    }
                    if(rparl < rmax && rperp < rmax){
                        xi2d[std::floor(rperp/binspacing)+nbins*std::floor(rparl/binspacing)].accumulate(value);
                    }
                }
            }
        }
        if(nfields > 10 && (ifield+1) % (nfields/10) == 0) {
            // Print status message in 10% intervals
            if(verbose) {
                std::cout << "Generating " << ifield+1 << "..." << std::endl;
            }
            // Save 1d stack to file every 10%
            if(snapshot) {
                std::string outFilename((boost::format("%s.snap%d.1d.dat") % prefix % int((ifield+1.)/nfields*10)).str());
                std::ofstream out(outFilename.c_str());
                boost::format outFormat("%.2f %.10f %.10f %d");
                for(int index = 0; index < nbins; ++index) {
                    out << (outFormat % ((index+.5)*binspacing)
                        % xi[index].mean() % xi[index].variance() % xi[index].count()) << std::endl;
                }
                out.close();
            }
        }
    }

    // Print extreme value mean, variance, and count
    if(verbose) {
        std::cout << boost::format("Extreme value mean, variance, count: %f %f %d")
            % extremeValues.mean() % extremeValues.variance() % extremeValues.count() << std::endl;
    }

    // Save 1d stack to file.
    {
        std::string outFilename(prefix+".1d.dat");
        std::ofstream out(outFilename.c_str());
        boost::format outFormat("%.2f %.10f %.10f %d");
        for(int index = 0; index < nbins; ++index) {
            out << (outFormat % ((index+.5)*binspacing)
                % xi[index].mean() % xi[index].variance() % xi[index].count()) << std::endl;
        }
        out.close();
    } 

    // Save 2d stack to file.
    {
        std::string outFilename(prefix+".2d.dat");
        std::ofstream out(outFilename.c_str());
        boost::format outFormat("%.2f %.2f %.10f %.10f %d");
        for(int index = 0; index < nbins*nbins; ++index) {
            out << (outFormat % ((index%nbins+.5)*binspacing) % ((index/nbins+.5)*binspacing)
                % xi2d[index].mean() % xi2d[index].variance() % xi2d[index].count()) << std::endl;
        }
        out.close();
    }

    // Done
    return 0;
}