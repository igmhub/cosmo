// Created 18-Feb-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/BroadbandPower.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "cosmo/RuntimeError.h"

#include "boost/bind.hpp"
#include "boost/ref.hpp"
#include "boost/foreach.hpp"

#include <cmath>

namespace local = cosmo;

local::BroadbandPower::BroadbandPower(double coef, double p,
double kmin, double rmin, double r0, double sigmaSq)
: _coef(coef), _p(p), _kmin(kmin), _rmin(rmin)
{
    if(kmin*rmin >= 1) throw RuntimeError("BroadbandPower: kmin*rmin < 1.");
    // Precompute constants
    double pi(4*std::atan(1));
    _twopi2 = 2*pi*pi;
    _kminp = std::pow(kmin,p);
    // Rescale coefficients to include "natural" normalization constants?
    if(r0 > 0) {
        if(sigmaSq <= 0) throw RuntimeError("BroadbandPower: expected sigmaSq > 0.");
        // Calculate the RMS fluctuations of PB(k,p)/B(p) within a sphere of radius r0
        PowerSpectrumPtr self(new PowerSpectrum(boost::bind(
            &BroadbandPower::operator(),boost::ref(*this),_1)));
        _coef = 1;
        double sigmaOld = getRmsAmplitude(self,r0);
        // Rescale coefs[n] so that a value of one would give the requested sigmaSq.
        double ratio(sigmaSq/(sigmaOld*sigmaOld));
        _coef = ratio*coef;
    }
}

local::BroadbandPower::~BroadbandPower() { }

double local::BroadbandPower::operator()(double kMpch) const {
    if(kMpch < 0) throw RuntimeError("BroadbandPower: expected wavenumber kMpch >= 0.");
    double krmin(kMpch*_rmin), k2(kMpch*kMpch);
    return _coef*(k2*kMpch/_twopi2)*std::exp(-krmin*krmin)/(_kminp+std::pow(kMpch,_p));
}

// explicit template instantiation for creating a function pointer to a TransferFunctionPowerSpectrum.

#include "likely/function_impl.h"

template local::PowerSpectrumPtr likely::createFunctionPtr<local::BroadbandPower>
    (boost::shared_ptr<local::BroadbandPower> pimpl);
