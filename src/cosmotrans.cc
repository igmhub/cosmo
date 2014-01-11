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
	double operator()(double k) const { return 1; }
};

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology multipole transforms");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
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

    cosmo::MultipoleTransform::Type type(cosmo::MultipoleTransform::SphericalBessel);
    boost::shared_ptr<Power> Pk(new Power());
    lk::GenericFunctionPtr PkPtr = lk::createFunctionPtr<Power>(Pk);

    int npts(100);
    std::vector<double> rvec(npts), xivec(npts);

    try {
    	cosmo::MultipoleTransform mt(type,0,0.1,10.,1e-3,40);
        mt.transform(PkPtr,rvec,xivec);
    }
    catch(std::runtime_error const &e) {
        std::cerr << "ERROR: exiting with an exception:\n  " << e.what() << std::endl;
        return -1;
    }

    return 0;
}