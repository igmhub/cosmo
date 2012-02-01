// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/RsdPowerSpectrumCorrelationFunction.h"

namespace local = cosmo;

local::RsdPowerSpectrumCorrelationFunction::RsdPowerSpectrumCorrelationFunction(
PowerSpectrumPtr powerSpectrumPtr, double rmin, double rmax, int nr) :
_monopole(powerSpectrumPtr, rmin, rmax, PowerSpectrumCorrelationFunction::Monopole, nr),
_quadrupole(powerSpectrumPtr, rmin, rmax, PowerSpectrumCorrelationFunction::Quadrupole, nr),
_hexadecapole(powerSpectrumPtr, rmin, rmax, PowerSpectrumCorrelationFunction::Hexadecapole, nr)
{
}

local::RsdPowerSpectrumCorrelationFunction::~RsdPowerSpectrumCorrelationFunction() { }

void local::RsdPowerSpectrumCorrelationFunction::setDistortion(double beta1, double beta2) {
    if(0 == beta2) beta2 = beta1;
    double betaSum(beta1+beta2), betaProd(beta1*beta2);
    _C0 = 1 + (1./3.)*betaSum + (1./5.)*betaProd;
    _C2 = (2./3.)*betaSum + (4./7.)*betaProd;
    _C4 = (8./35.)*betaProd;
}

double local::RsdPowerSpectrumCorrelationFunction::operator()(double rMpch, double mu) const {
    // Evaluate the l=2,4 Legendre polynomials in mu
    double mu2(mu*mu);
    double P2(1.5*mu2-0.5), P4(4.375*mu2*mu2-3.75*mu2+0.375);
    // Combine the l=0,2,4 pieces.
    return _C0*_monopole(rMpch) + P2*_C2*_quadrupole(rMpch) + P4*_C4*_hexadecapole(rMpch);
}