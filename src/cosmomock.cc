// Created 06-Aug-2013 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
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

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    std::string rvectors,kvectors,outfile;
    po::options_description cli("Mock generator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("rvectors", po::value<std::string>(&rvectors)->default_value(""),
            "Filename to read r-vectors from")
        ("kvectors", po::value<std::string>(&kvectors)->default_value(""),
            "Filename to read k-vectors from")
        ("output,o", po::value<std::string>(&outfile)->default_value("mock.dat"),
            "Filename to save generated mock to")
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
    bool verbose(vm.count("verbose")),rmu(vm.count("rmu"));

    // Read the r-vectors file
    if(0 == rvectors.length()) {
        std::cerr << "Missing rvectors parameter." << std::endl;
        return -2;
    }
    std::vector<std::vector<double> > rvec(3);
    try {
        std::ifstream in(rvectors.c_str());
        lk::readVectors(in,rvec);
        in.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while reading " << rvectors << ": " << e.what() << std::endl;
        return -3;
    }
    int npixels = rvec[0].size();
    if(verbose) {
        std::cout << "Read " << npixels << " k-vectors from " << rvectors
            << std::endl;
    }

    // Read the k-vectors file
    if(0 == kvectors.length()) {
        std::cerr << "Missing kvectors parameter." << std::endl;
        return -2;
    }
    std::vector<std::vector<double> > kvec(5);
    try {
        std::ifstream in(kvectors.c_str());
        lk::readVectors(in,kvec);
        in.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while reading " << kvectors << ": " << e.what() << std::endl;
        return -3;
    }
    int nmodes = kvec[0].size();
    if(verbose) {
        std::cout << "Read " << nmodes << " k-vectors from " << kvectors
            << std::endl;
    }

    std::ofstream out(outfile.c_str());

    // Evaluate realization (kvec) at each survey pixel (rvec)
    double wgt = 1;
    for(int i = 0; i < npixels; ++i) {
        double delta(0);
        for(int j = 0; j < nmodes; ++j) {
            double dot = kvec[0][j]*rvec[0][i]+kvec[1][j]*rvec[1][i]+kvec[2][j]*rvec[2][i];
            delta += 2*kvec[3][j]*std::cos(dot+kvec[4][j]);
        }
        out << rvec[0][i] << ' ' << rvec[1][i] << ' ' << rvec[2][i] << ' ' << delta << ' ' << wgt << std::endl;
    }

    out.close();

    return 0;
}
