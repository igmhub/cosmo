// Created 10-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/OneDimensionalPowerSpectrum.h"
#include "cosmo/RuntimeError.h"

#include "likely/Integrator.h"
#include "likely/Interpolator.h"

#include "boost/bind.hpp"
#include "boost/math/special_functions/bessel.hpp"

#include <cmath>

namespace local = cosmo;

local::OneDimensionalPowerSpectrum::OneDimensionalPowerSpectrum(
PowerSpectrumPtr powerSpectrum, double radius, double kmin, double kmax, int nk)
: _powerSpectrum(powerSpectrum), _radius(radius), _kmin(kmin), _kmax(kmax), _nk(nk)
{
    if(kmin <= 0) {
        throw RuntimeError("OneDimensionalPowerSpectrum: invalid kmin <= 0.");
    }
    if(kmax <= kmin) {
        throw RuntimeError("OneDimensionalPowerSpectrum: invalid kmax <= kmin.");
    }
    if(nk < 2) {
        throw RuntimeError("OneDimensionalPowerSpectrum: invalid nk < 2.");
    }    
}

local::OneDimensionalPowerSpectrum::~OneDimensionalPowerSpectrum() { }

double local::OneDimensionalPowerSpectrum::operator()(double kMpch) const {
    if(kMpch < _kmin) {
        throw RuntimeError("PowerSpectrumCorrelationFunction: k < kmin.");
    }
    if(kMpch > _kmax) {
        throw RuntimeError("PowerSpectrumCorrelationFunction: k > kmax.");
    }
    if(!_interpolator) {
        // Create an integrator.
        likely::Integrator::IntegrandPtr integrand(new likely::Integrator::Integrand(
            boost::bind(&OneDimensionalPowerSpectrum::_integrand,this,_1)));
        likely::Integrator integrator(integrand,1e-7,1e-6);
        // Allocate temporary space for the interpolation tables.
        likely::Interpolator::CoordinateValues logkValues(_nk), pValues(_nk);
        // Loop over logarithmic steps in r to build the interpolation tables.
        double logkmin(std::log(_kmin)), logkmax(std::log(_kmax));
        double dlogk((logkmax-logkmin)/(_nk-1));
        logkValues[_nk-1] = logkmax;
        double kzLast(_kmax);
        _kz2 = _kmax*_kmax;
        pValues[_nk-1] = integrator.integrateUp(_kmax);
        /** Introduce a Jeans-scale cutoff in the r=0 limit *
        if(_radius == 0) {
            double pi(4*atan(1));
            std::cout << "kcut = " << pi/0.1 << std::endl;
            pValues[_nk-1] -= integrator.integrateUp(pi/0.1);
        }
        **/
        for(int i = _nk-2; i >= 0; --i) {
            double logk = logkmin + i*dlogk;
            logkValues[i] = logk;
            double kz(std::exp(logk));
            if(_radius == 0) {
                pValues[i] = pValues[i+1] + integrator.integrateSmooth(kz,kzLast);
                kzLast = kz;
            }
            else {
                _kz2 = kz*kz;
                pValues[i] = integrator.integrateUp(kz);
            }
        }
        for(int i = 0; i < _nk; ++i) {
            double kz(std::exp(logkValues[i]));
            pValues[i] *= kz;
        }
        _interpolator.reset(new likely::Interpolator(logkValues,pValues,"cspline"));
    }
    return (*_interpolator)(std::log(kMpch));
}

double local::OneDimensionalPowerSpectrum::_integrand(double kval) const {
    double k2(kval*kval);
    double result((*_powerSpectrum)(kval)/k2);
    if(_radius > 0) {
        double ktr(std::sqrt(k2-_kz2)*_radius);
        double wfun(2*boost::math::cyl_bessel_j(1,ktr)/ktr);
        result *= wfun*wfun;
    }
    else if(_radius < 0) {
        double ktr2((k2-_kz2)*_radius*_radius);
        double wfun(std::exp(-ktr2/2));
        result *= wfun*wfun;
    }
    return result;
}
