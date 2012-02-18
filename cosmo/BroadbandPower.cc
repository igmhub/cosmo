// Created 18-Feb-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/BroadbandPower.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "cosmo/RuntimeError.h"

#include "boost/bind.hpp"
#include "boost/ref.hpp"
#include "boost/foreach.hpp"

#include <cmath>

//!!
#include <iostream>

namespace local = cosmo;

local::BroadbandPower::BroadbandPower(int nmin, std::vector<double> coefs,
double rmin, double rmax, double r0, double sigmaSq)
: _nmin(nmin), _nmax(nmin+coefs.size()), _coefs(coefs), _rmin(rmin), _rmax(rmax)
{
    if(rmax <= rmin) throw RuntimeError("BroadbandPower: expected rmax > rmin.");
    // Rescale coefficients to include "natural" normalization constants?
    if(r0 > 0) {
        if(sigmaSq <= 0) throw RuntimeError("BroadbandPower: expected sigmaSq > 0.");
        for(int dn = 0; dn < coefs.size(); ++dn) {
            // Calculate the RMS fluctuations of PB(k,n) within a sphere of radius r0
            PowerSpectrumPtr PBptr(new PowerSpectrum(boost::bind(
                &BroadbandPower::evaluatePB,boost::ref(*this),_1,_nmin+dn)));
            double sigmaOld = getRmsAmplitude(PBptr, r0);
            // Rescale coefs[n] so that a value of one would give the requested sigmaSq.
            double ratio(sigmaSq/(sigmaOld*sigmaOld));
            std::cout << _nmin + dn << " => " << ratio << std::endl;
            _coefs[dn] *= ratio;            
        }
    }
    for(int dn = 0; dn < coefs.size(); ++dn) _powrmax.push_back(std::pow(rmax,nmin+dn));
}

local::BroadbandPower::~BroadbandPower() { }

double local::BroadbandPower::evaluatePB(double k, int n) const {
    if(n < _nmin || n >= _nmax) throw RuntimeError("BroadbandPower: exponent n out of range.");
    double krmin(k*_rmin), krmin2(krmin*krmin), krmax(k*_rmax);
    return std::exp(-krmin2)*_powrmax[n-_nmin]/(1+std::pow(krmax,n));
}

double local::BroadbandPower::operator()(double kMpch) const {
    double result(0);
    for(int dn = 0; dn < _coefs.size(); ++dn) {
        result += _coefs[dn]*evaluatePB(kMpch,_nmin+dn);
    }
    return result;
}