// Created 05-Aug-2013 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Correlation function estimator.

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

namespace po = boost::program_options;
namespace lk = likely;

void bruteForce(std::vector<std::vector<double> > const &columns) {
    int n(columns[0].size()), npair(0), nused(0);
    for(int i = 0; i < n-1; ++i) {
        double xi = columns[0][i];
        double yi = columns[1][i];
        double zi = columns[2][i];
        for(int j = i+1; j < n; ++j) {
            double dx = xi - columns[0][j];
            double dy = yi - columns[1][j];
            double dz = zi - columns[2][j];
            double r = std::sqrt(dx*dx+dy*dy+dz*dz);
            npair++;
            if(r < 200) nused++;
        }
    }
    std::cout << "used " << nused << " of " << npair << " pairs." << std::endl;
}

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    std::string infile;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read field samples from")
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

    // Read the input file
    if(0 == infile.length()) {
        std::cerr << "Missing infile parameter." << std::endl;
        return -2;
    }
    std::vector<std::vector<double> > columns(5);
    try {
        std::ifstream in(infile.c_str());
        lk::readVectors(in,columns);
        in.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while reading " << infile << ": " << e.what() << std::endl;
        return -3;
    }
    if(verbose) {
        std::cout << "Read " << columns[0].size() << " rows from " << infile
            << std::endl;
    }

    bruteForce(columns);

    return 0;
}
