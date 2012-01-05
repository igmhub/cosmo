// Created 04-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Converts an (ra,dec,z) catalog into an (x,y,z) catalog in Mpc/h.

#include "cosmo/cosmo.h"

#include "boost/program_options.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology calculator");
    double OmegaLambda,OmegaMatter,OmegaBaryon,hubbleConstant,cmbTemp,spectralIndex,
        zval,kval,kmin,kmax,rval,rmin,rmax;
    int nk,nr;
    std::string inputName, outputName;
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.728),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0),
            "Present-day value of OmegaMatter or zero for 1-OmegaLambda.")
        ("input-name,i", po::value<std::string>(&inputName)->default_value(""),
            "Name of the input file containing ra,dec,z values to read.")
        ("output-name,o", po::value<std::string>(&outputName)->default_value(""),
            "Name of the output file containing x,y,z values to write.")
        ("bounds", "Calculates and prints bounding box of converted points.")
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
    bool verbose(vm.count("verbose")), bounds(vm.count("bounds"));

    if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
    cosmo::AbsHomogeneousUniversePtr cosmology(
        new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));
    
    std::ofstream out(outputName.c_str());
    std::ifstream in(inputName.c_str());
    double ra,dec,z;
    int count(0);
    double deg2rad(atan2(1,0)/90);
    double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
    while(in.good() && !in.eof()) {
        in >> ra >> dec >> z;
        if(!in.good() || in.eof()) break;
        if(z <= 0) {
            std::cerr << "Bad redshift on line " << (count+1) << " : z = " << z << std::endl;
        }
        else {
            double s(cosmology->getLineOfSightComovingDistance(z));
            double RA(deg2rad*ra), DEC(deg2rad*dec);
            double cosDEC(std::cos(DEC));
            double X(s*cosDEC*std::cos(RA)), Y(s*cosDEC*std::sin(RA)), Z(s*std::sin(DEC));
            out << X << ' ' << Y << ' ' << Z << std::endl;
            if(bounds) {
                if(0 == count) {
                    Xmin = Xmax = X;
                    Ymin = Ymax = Y;
                    Zmin = Zmax = Z;
                }
                else {
                    if(X < Xmin) Xmin = X;
                    else if(X > Xmax) Xmax = X;
                    if(Y < Ymin) Ymin = Y;
                    else if(Y > Ymax) Ymax = Y;
                    if(Z < Zmin) Zmin = Z;
                    else if(Z > Zmax) Zmax = Z;
                }
            }
        }
        count++;
    }
    if(bounds) {
        std::cout << "Bounding box [" << Xmin << ',' << Xmax << "] x [" << Ymin << ',' << Ymax
            << "] x [" << Zmin << ',' << Zmax << "]" << std::endl;
    }
    if(verbose) {
        std::cout << "Converted " << count << " lines." << std::endl;
    }
    in.close();
    out.close();
    return 0;
}