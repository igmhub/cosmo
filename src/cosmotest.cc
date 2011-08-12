// Created 10-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/cosmo.h"

#include "boost/math/special_functions/pow.hpp"
#include "boost/math/special_functions/expint.hpp"

#include <iostream>
#include <cmath>

// Calculates k^3/(2pi^2) P(k) = k for an input wavenumber in 1/(Mpc/h).
double powerSpectrum(double kval) {
    return kval;
}

double powerSpectrum2(double kval) {
    return 0.050660591821168885722; // 1/(2pi^2)
}

int main(int argc, char **argv) {

    double pi(4*std::atan(1)), rootpi(std::sqrt(pi));

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
    
    // Calculate the 1D power spectrum from kz = 0.001 - 10 /(Mpc/h)
    double kmin(0.001), kmax(10), radius(1);
    double r3(boost::math::pow<3>(radius));
    cosmo::OneDimensionalPowerSpectrum oned(powerPtr,-radius,kmin,kmax);
    for(int i = 0; i < ntest; ++i) {
        double kz(kmin*std::pow(kmax/kmin,i/(ntest-1.)));
        double pval(oned(kz));
        double kzr(kz*radius);
        double exact(-kz/2*std::exp(kzr*kzr)*boost::math::expint(-kzr*kzr));
        std::cout << "P1D(kz = " << kz << " /(Mpc/h)) = " << pval << " (error = "
            << pval-exact << ")" << std::endl;
    }
    
    // Check the normalization of JMLG's paper draft in his Eqn.1 and Fig.1
    cosmo::PowerSpectrumPtr powerPtr2(new cosmo::PowerSpectrum(powerSpectrum2));
    cosmo::OneDimensionalPowerSpectrum oned2(powerPtr2,0,kmin,kmax);
    double kval(1);
    std::cout << (pi/kval)*oned2(kval) << std::endl;
    
    return 0;
}
