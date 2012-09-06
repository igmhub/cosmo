// Created 09-Jul-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Generates a Gaussian random field.

// $ time ./cosmogrf --spacing 5 --load-power ../../baofit/models/DR9LyaMocks_matterpower.dat --nx 512 
//      --powerfile power.dat --npairs 1000000000 --corrfile corr.dat
// real    7m29.486s
// user    7m28.809s
// sys 0m0.680s

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <iostream>
#include <fstream>
#include <string>

namespace po = boost::program_options;
namespace lk = likely;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    double spacing;
    long npairs;
    int nx,ny,nz,seed,pairseed,nbins,nkbins;
    std::string loadPowerFile, corrfile, powerfile;
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
            "Random seed to use for GRF.")
        ("corrfile", po::value<std::string>(&corrfile)->default_value(""),
            "Name of correlation function output file, leave blank to skip.")
        ("npairs", po::value<long>(&npairs)->default_value(1000000),
            "Number of pairs to use in correlation function estimate.")
        ("pair-seed", po::value<int>(&pairseed)->default_value(42),
            "Random seed to use for correlation function pairs.")
        ("nbins", po::value<int>(&nbins)->default_value(50),
            "Number of r bins to use for correlation function measurement.")
        ("powerfile", po::value<std::string>(&powerfile)->default_value(""),
            "Name of power spectrum output file, leave blank to skip.")
        ("nkbins", po::value<int>(&nkbins)->default_value(100),
            "Number of k bins to use for power spectrum measurement.")
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
    if(nbins <= 0) {
        std::cerr << "nbins must be > 0" << std::endl;
        return -2;
    }    
    if(nkbins <= 0) {
        std::cerr << "nkbins must be > 0" << std::endl;
        return -3;
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

    // Generate delta field in k-space.
    generator.generateFieldK();

    // Perform power spectrum estimate from k-space delta field.
    if (powerfile.length() > 0) {
        // Prepare power spectrum bins.
        double kmax = pi/spacing, kmin = pi/(spacing*std::pow(nx*ny*nz,1./3.));
        double binsize = (kmax-kmin)/nkbins;
        boost::scoped_array<lk::WeightedAccumulator> 
            powerAccumulator(new lk::WeightedAccumulator[nkbins]);
        // k-space grid spacing.
        double dkx(2*pi/(nx*spacing)), dky(2*pi/(ny*spacing)), dkz(2*pi/(nz*spacing));
        double powerNorm = nx*ny*nz*spacing*spacing*spacing;
        for(int ix = 0; ix < nx; ++ix){
            double kx = (ix > nx/2 ? ix-nx : ix)*dkx;
            for(int iy = 0; iy < ny; ++iy){
                double ky = (iy > ny/2 ? iy-ny : iy)*dky;
                for(int iz = 0; iz < nz; ++iz){
                    double kz = (iz > nz/2 ? iz-nz : iz)*dkz;
                    double ksq = kx*kx + ky*ky + kz*kz;
                    double k = std::sqrt(ksq);
                    int index = (int) ((k-kmin)/binsize);
                    if (index < 0 || index >= nkbins) continue;
                    double redk(generator.getFieldKRe(ix,iy,iz)), imdk(generator.getFieldKIm(ix,iy,iz));
                    powerAccumulator[index].accumulate(powerNorm*(redk*redk+imdk*imdk));
                }
            }
        }
        // Output power spectrum.
        std::ofstream out(powerfile.c_str());
        for(int i = 0; i < nkbins; ++i) {
            out << boost::format("%d %f %f %f %d") 
                % i % (kmin + (i+.5)*binsize)
                % powerAccumulator[i].mean() % powerAccumulator[i].variance() 
                % powerAccumulator[i].count() << std::endl;
        }
        out.close();
    }

    // Perform FFT to realspace.
    generator.transformFieldToR();
    
    // Perform power spectrum estimate from k-space delta field.
    if (corrfile.length() > 0) {
        double rmin(0), rmax(200);
        double binsize = (rmax - rmin)/nbins;
        // Prepare correlation function accumulators.
        boost::scoped_array<lk::WeightedAccumulator> 
            corrdidj(new lk::WeightedAccumulator[nbins]),
            corrdi(new lk::WeightedAccumulator[nbins]),
            corrdj(new lk::WeightedAccumulator[nbins]);
        // Initialize the random number source for pair generating.
        lk::RandomPtr random = lk::Random::instance();
        random->setSeed(pairseed);
        for(long ipair = 0; ipair < npairs; ++ipair) {
            // First random position.
            long i = random->getInteger(1,nx*ny*nz);
            int ix(i % nx), iy(i/nx % ny), iz(i/(nx*ny) % nz);
            // Second random position chosen from to be within a sphere 
            // of radius rmax from the first position.
            double jr = rmax/spacing*std::pow(random->getUniform(),1/3.);
            double jphi = 2*pi*random->getUniform();
            double jtheta = std::acos(2*random->getUniform()-1);
            int jx = ix + (int) jr*cos(jphi)*sin(jtheta);
            int jy = iy + (int) jr*sin(jphi)*sin(jtheta);
            int jz = iz + (int) jr*cos(jtheta);
            if( jx < 0 || jx >= nx || jy < 0 || jy >= ny || jz < 0 || jz >= nz ) {
                ipair--;
                continue;
            }
            // Calculate distance between positions.
            double dx(ix-jx), dy(iy-jy), dz(iz-jz);
            double r = spacing*std::sqrt(dx*dx+dy*dy+dz*dz);
            int index = (int) ((r-rmin)/binsize);
            if(index >= nbins || index < 0 ) {
                ipair--;
                continue;
            }
            // Accumulate values.
            double idelta = generator.getField(ix,iy,iz);
            double jdelta = generator.getField(jx,jy,jz);
            corrdidj[index].accumulate(idelta*jdelta);
            corrdi[index].accumulate(idelta);
            corrdj[index].accumulate(jdelta);
        }
        // Save correlation function to file.
        std::ofstream out(corrfile.c_str());
        for(int i = 0; i < nbins; ++i) {
            out << boost::format("%d %f %f %f %f %f %f %f %d") 
                % i % (rmin + (i+.5)*binsize)
                % corrdidj[i].mean() % corrdidj[i].variance() 
                % corrdi[i].mean() % corrdi[i].variance() 
                % corrdj[i].mean() % corrdj[i].variance() 
                % corrdidj[i].count()  << std::endl;
        }
        out.close();
    }

    if(verbose) {
        // Calculate the statistics of the generated delta field.
        lk::WeightedAccumulator accumulator;
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
    }
    
    return 0;
}
