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
PowerSpectrumPtr powerSpectrum, double rmin, double rmax, Multipole multipole, int nr)
: _powerSpectrum(powerSpectrum), _rmin(rmin), _rmax(rmax), _multipole(multipole), _nr(nr)
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
        // Create separate integrators for k <= 4*pi/r and k > 4*pi/r.
        likely::Integrator::IntegrandPtr integrand1(new likely::Integrator::Integrand(
            boost::bind(&PowerSpectrumCorrelationFunction::_integrand1,this,_1)));
        likely::Integrator integrator1(integrand1,1e-8,1e-6);        
        likely::Integrator::IntegrandPtr integrand2(new likely::Integrator::Integrand(
            boost::bind(&PowerSpectrumCorrelationFunction::_integrand2,this,_1)));
        likely::Integrator integrator2(integrand2,1e-7,0); // sin-integrand needs more accuracy than cos      
        likely::Integrator::IntegrandPtr integrand3(new likely::Integrator::Integrand(
            boost::bind(&PowerSpectrumCorrelationFunction::_integrand3,this,_1)));
        likely::Integrator integrator3(integrand3,1e-6,0);        
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
            double kcut(4*pi/_radius);
            xiValues[i] = integrator1.integrateSingular(0,kcut) +
                integrator2.integrateOscUp(kcut,_radius,true);
            // Add the cosine oscillating part for l = 2,4
            if(_multipole != Monopole) xiValues[i] += integrator3.integrateOscUp(kcut,_radius,false);
        }
        _interpolator.reset(new likely::Interpolator(logrValues,xiValues,"cspline"));
    }
    return (*_interpolator)(std::log(rMpch));
}

double local::PowerSpectrumCorrelationFunction::_integrand1(double kval) const {
    double kr(kval*_radius),kr2(kr*kr);
    double result = (*_powerSpectrum)(kval);
    switch(_multipole) {
    case Monopole:
        result *= boost::math::sinc_pi(kr)/kval;
        break;
    case Quadrupole:
        result *= (std::sin(kr)*(kr2 - 3) + std::cos(kr)*(3*kr))/(kr2*kr*kval);
        break;
    case Hexadecapole:
        result *= (std::sin(kr)*(kr2*kr2 - 45*kr2 + 105) +
            std::cos(kr)*(10*kr2 - 105)*kr)/(kr2*kr2*kr*kval);
        break;
    }
    return result;
}

double local::PowerSpectrumCorrelationFunction::_integrand2(double kval) const {
    double kr(kval*_radius),kr2(kr*kr);
    double result = (*_powerSpectrum)(kval);
    switch(_multipole) {
    case Monopole:
        result *= 1/(kr*kval);
        break;
    case Quadrupole:
        result *= (kr2 - 3)/(kr2*kr*kval);
        break;
    case Hexadecapole:
        result *= (kr2*kr2 - 45*kr2 + 105)/(kr2*kr2*kr*kval);
        break;
    }
    return result;
}

double local::PowerSpectrumCorrelationFunction::_integrand3(double kval) const {
    double kr(kval*_radius),kr2(kr*kr);
    double result = (*_powerSpectrum)(kval);
    switch(_multipole) {
    case Monopole:
        result = 0;
        break;
    case Quadrupole:
        result *= 3/(kr2*kval);
        break;
    case Hexadecapole:
        result *= (10*kr2 - 105)/(kr2*kr2*kval);
        break;
    }
    return result;
}
