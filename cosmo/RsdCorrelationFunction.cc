// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/RsdCorrelationFunction.h"

namespace local = cosmo;

local::RsdCorrelationFunction::RsdCorrelationFunction(
CorrelationFunctionPtr xi0, CorrelationFunctionPtr xi2, CorrelationFunctionPtr xi4) :
_xi0(xi0), _xi2(xi2), _xi4(xi4)
{
}

local::RsdCorrelationFunction::~RsdCorrelationFunction() { }

void local::RsdCorrelationFunction::setDistortion(double beta1, double beta2) {
    if(0 == beta2) beta2 = beta1;
    double betaSum(beta1+beta2), betaProd(beta1*beta2);
    _C0 = 1 + (1./3.)*betaSum + (1./5.)*betaProd;
    _C2 = (2./3.)*betaSum + (4./7.)*betaProd;
    _C4 = (8./35.)*betaProd;
}

double local::RsdCorrelationFunction::operator()(double rMpch, double mu) const {
    // Evaluate the l=2,4 Legendre polynomials in mu
    double mu2(mu*mu);
    double P2(1.5*mu2-0.5), P4(4.375*mu2*mu2-3.75*mu2+0.375);
    // Combine the l=0,2,4 pieces.
    return _C0*(*_xi0)(rMpch) + P2*_C2*(*_xi2)(rMpch) + P4*_C4*(*_xi4)(rMpch);
}

double local::RsdCorrelationFunction::operator()(double rMpch, Multipole multipole) const {
    switch(multipole) {
    case Hexadecapole:
        return (*_xi4)(rMpch);
    case Quadrupole:
        return (*_xi2)(rMpch);
    default:
        return (*_xi0)(rMpch);
    }
}
