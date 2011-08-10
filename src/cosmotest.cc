// Created 10-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/cosmo.h"

#include <iostream>
#include <cmath>

// Calculates k^3/(2pi^2) P(k) for an input wavenumber in 1/(Mpc/h)
double powerSpectrum(double kval) {
    return kval;
}

int main(int argc, char **argv) {

    double pi(4*std::atan(1));

    // Create a shared pointer to a power spectrum function object.
    cosmo::PowerSpectrumPtr powerPtr(new cosmo::PowerSpectrum(powerSpectrum));
    
    // Calculate sigma8 with a top-hat weight function.
    double r8(8);
    double exactSigma8 = std::sqrt(3*pi/(5*r8));
    double sigma8 = cosmo::getRmsAmplitude(powerPtr,r8);
    std::cout << "sigma(8 Mpc/h) = " << sigma8 << " (error = " << sigma8-exactSigma8
        << ")"<< std::endl;

    // Calculate sigmaQSO with r = 0.04 Mpc/h and a Gaussian weighting function.
    double rQSO(0.04);
    double exactSigmaQSO = std::pow(pi,1./4.)/std::sqrt(2*rQSO);
    double sigmaQSO = cosmo::getRmsAmplitude(powerPtr,rQSO,true);
    std::cout << "sigmaQSO = " << sigmaQSO << " (error = " << sigmaQSO-exactSigmaQSO
        << ")" << std::endl;
        
    // Calculate the correlation function from r = 0.01 - 100 Mpc/h.
    double rmin(0.01), rmax(100);
    cosmo::PowerSpectrumCorrelationFunction xi(powerPtr,rmin,rmax);
    int ntest(20);
    for(int i = 0; i < ntest; ++i) {
        double rval(rmin*std::pow(rmax/rmin,i/(ntest-1.)));
        double xival(xi(rval));
        double exact(pi/(2*rval));
        std::cout << "xi(r = " << rval << " Mpc/h) = " << xival << " (error = "
            << xival-exact << ")" << std::endl;
    }
    return 0;
}
