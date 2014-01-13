// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/cosmo.h"
#include "likely/likely.h"
#include "likely/function_impl.h"

#include "boost/program_options.hpp"
#include "boost/bind.hpp"

namespace po = boost::program_options;
namespace lk = likely;

class Power {
public:
	Power() { }
	double operator()(double k) const { return k; }
};

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology multipole transforms");
    int ell,minSamplesPerCycle,minSamplesPerDecade;
    double min,max,veps;
    cli.add_options()
        ("help,h", "prints this info and exits.")
        ("verbose", "prints additional information.")
        ("hankel", "performs a Hankel transform (default is spherical Bessel")
        ("ell", po::value<int>(&ell)->default_value(0),
            "multipole number of transform to calculate")
        ("min", po::value<double>(&min)->default_value(0.1),
            "minimum value of transformed coordinate")
        ("max", po::value<double>(&max)->default_value(10.),
            "maximum value of transformed coordinate")
        ("veps", po::value<double>(&veps)->default_value(1e-3),
            "desired transform accuracy")
        ("measure", "does initial measurements to optimize FFT plan")
        ("min-samples-per-cycle", po::value<int>(&minSamplesPerCycle)->default_value(2),
            "minimum number of samples per cycle to use for transform convolution")
        ("min-samples-per-decade", po::value<int>(&minSamplesPerDecade)->default_value(40),
            "minimum number of samples per decade to use for transform convolution")
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

    cosmo::MultipoleTransform::Type ttype(hankel ?
        cosmo::MultipoleTransform::Hankel :
        cosmo::MultipoleTransform::SphericalBessel);

    cosmo::MultipoleTransform::Strategy strategy(measure ?
        cosmo::MultipoleTransform::MeasurePlan :
        cosmo::MultipoleTransform::EstimatePlan);

    boost::shared_ptr<Power> Pk(new Power());
    lk::GenericFunctionPtr PkPtr = lk::createFunctionPtr<Power>(Pk);

    int npts(100);
    std::vector<double> rvec(npts), xivec(npts);

    try {
    	cosmo::MultipoleTransform mt(ttype,ell,min,max,veps,strategy,
            minSamplesPerCycle,minSamplesPerDecade);
        std::vector<double> const& ugrid = mt.getUGrid(), vgrid = mt.getVGrid();
        if(verbose) {
            std::cout << "Truncation fraction is " << mt.getTruncationFraction() << std::endl;
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
        for(int i = 0; i < results.size(); ++i) {
            std::cout << results[i] << ", ";
        }
        std::cout << std::endl;
    }
    catch(std::runtime_error const &e) {
        std::cerr << "ERROR: exiting with an exception:\n  " << e.what() << std::endl;
        return -1;
    }

    return 0;
}