// Created 05-Aug-2013 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Correlation function estimator.

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
    po::options_description cli("Correlation function estimator");
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
    
    return 0;
}
