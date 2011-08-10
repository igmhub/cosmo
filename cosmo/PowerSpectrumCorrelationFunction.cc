// Created 09-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/PowerSpectrumCorrelationFunction.h"
#include "cosmo/RuntimeError.h"

#include "likely/Integrator.h"
#include "likely/Interpolator.h"

#include "boost/bind.hpp"
#include "boost/math/special_functions/sinc.hpp"

#include <cmath>

namespace local = cosmo;

local::PowerSpectrumCorrelationFunction::PowerSpectrumCorrelationFunction(
PowerSpectrumPtr powerSpectrum, double rmin, double rmax, int nr)
: _powerSpectrum(powerSpectrum), _rmin(rmin), _rmax(rmax), _nr(nr)
{
    if(rmin <= 0) {
        throw RuntimeError("PowerSpectrumCorrelationFunction: invalid rmin <= 0.");
    }
    if(rmax <= rmin) {
        throw RuntimeError("PowerSpectrumCorrelationFunction: invalid rmax <= rmin.");
    }
    if(nr < 2) {
        throw RuntimeError("PowerSpectrumCorrelationFunction: invalid nr < 2.");
    }
}

local::PowerSpectrumCorrelationFunction::~PowerSpectrumCorrelationFunction() { }

double local::PowerSpectrumCorrelationFunction::operator()(double rMpch) const {
    if(rMpch < _rmin) {
        throw RuntimeError("PowerSpectrumCorrelationFunction: r < rmin.");
    }
    if(rMpch > _rmax) {
        throw RuntimeError("PowerSpectrumCorrelationFunction: r > rmax.");
    }
    if(!_interpolator) {
        // Create separate integrators for k <= pi/r and k > pi/r.
        likely::Integrator::IntegrandPtr integrand1(new likely::Integrator::Integrand(
            boost::bind(&PowerSpectrumCorrelationFunction::_integrand1,this,_1)));
        likely::Integrator integrator1(integrand1,1e-8,1e-6);        
        likely::Integrator::IntegrandPtr integrand2(new likely::Integrator::Integrand(
            boost::bind(&PowerSpectrumCorrelationFunction::_integrand2,this,_1)));
        likely::Integrator integrator2(integrand2,1e-6,0);        
        // Allocate temporary space for the interpolation tables.
        likely::Interpolator::CoordinateValues logrValues(_nr), xiValues(_nr);
        // Loop over logarithmic steps in r to build the interpolation tables.
        double logrmin(std::log(_rmin)), logrmax(std::log(_rmax));
        double dlogr((logrmax-logrmin)/(_nr-1));
        double pi(4*std::atan(1));
        for(int i = 0; i < _nr; ++i) {
            double logr = logrmin + i*dlogr;
            logrValues[i] = logr;
            // Save the current radius so it can be accessed by our integrand methods.
            _radius = std::exp(logr);
            // Calculate the integral over each domain separately.
            double kcut(pi/_radius);
            xiValues[i] = integrator1.integrateSingular(0,kcut) +
                integrator2.integrateOscUp(kcut,_radius);
        }
        _interpolator.reset(new likely::Interpolator(logrValues,xiValues,"cspline"));
    }
    return (*_interpolator)(std::log(rMpch));
}

double local::PowerSpectrumCorrelationFunction::_integrand1(double kval) const {
    double kr(kval*_radius);
    return (*_powerSpectrum)(kval)*boost::math::sinc_pi(kr)/kval;
}

double local::PowerSpectrumCorrelationFunction::_integrand2(double kval) const {
    double kr(kval*_radius);
    return (*_powerSpectrum)(kval)/(kr*kval);
}
