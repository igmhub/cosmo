// Created 20-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/DistortedPowerCorrelation.h"
#include "cosmo/AdaptiveMultipoleTransform.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "cosmo/RuntimeError.h"

#include "likely/Interpolator.h"

#include "boost/foreach.hpp"
#include "boost/bind.hpp"

namespace local = cosmo;

local::DistortedPowerCorrelation::DistortedPowerCorrelation(likely::GenericFunctionPtr power,
RMuFunctionCPtr distortion, double rmin, double rmax, int nr,
int ellMax, bool symmetric, double relerr, double abserr, double abspow)
: _power(power), _distortion(distortion), _ellMax(ellMax), _symmetric(symmetric),
_relerr(relerr), _abserr(abserr), _abspow(abspow), _initialized(false)
{
	if(rmax <= rmin) {
		throw RuntimeError("DistortedPowerCorrelation: expected rmin < rmax.");
	}
	if(rmin <= 0) {
		throw RuntimeError("DistortedPowerCorrelation: expected rmin > 0.");
	}
	if(nr < 2) {
		throw RuntimeError("DistortedPowerCorrelation: expected nr >= 2.");
	}
	if(ellMax < 0) {
		throw RuntimeError("DistortedPowerCorrelation: expected ellMax >= 0.");
	}
	if(symmetric && (ellMax%2 == 1)) {
		throw RuntimeError("DistortedPowerCorrelation: expected even ellMax when symmetric.");
	}
	// initialize the r grid we will use for interpolation
	_rgrid.reserve(nr);
	double dr = (rmax - rmin)/(nr-1.);
	for(int i = 0; i < nr; ++i) {
		_rgrid.push_back(rmin + dr*i);
	}
	// create a transform object for each moment
	int dell = symmetric ? 2 : 1;
	int nell = 1+ellMax/dell;
	_transformer.reserve(nell);
	_xiMoments.reserve(nell);
	_interpolator.reserve(nell);
	for(int ell = 0; ell <= ellMax; ell += dell) {
		double coef = multipoleTransformNormalization(ell,3,+1);
		// use the same relerr for each ell and share abserr equally
		AdaptiveMultipoleTransformPtr amt(new AdaptiveMultipoleTransform(
			MultipoleTransform::SphericalBessel,ell,coef,_rgrid,relerr,abserr/nell,abspow));
		_transformer.push_back(amt);
		_xiMoments.push_back(std::vector<double>(nr,0.));
		_interpolator.push_back(likely::InterpolatorPtr());
	}
}

local::DistortedPowerCorrelation::~DistortedPowerCorrelation() { }

double local::DistortedPowerCorrelation::getPower(double k, double mu) const {
	if(mu < 0 || mu > 1) {
		throw RuntimeError("DistortedPowerCorrelation::getPower: expected 0 <= mu <= 1.");
	}
	return (*_power)(k)*(*_distortion)(k,mu);
}

double local::DistortedPowerCorrelation::getPowerMultipole(double k, int ell) const {
	if(ell < 0 || ell > _ellMax || (_symmetric && (ell%2))) {
		throw RuntimeError("DistortedPowerCorrelation::getPowerMultipole: invalid ell.");
	}
	likely::GenericFunctionPtr fOfMuPtr(
		new likely::GenericFunction(boost::bind(
			&DistortedPowerCorrelation::getPower,this,k,_1)));
	return getMultipole(fOfMuPtr, ell);
}

double local::DistortedPowerCorrelation::getCorrelationMultipole(double r, int ell) const {
	if(!_initialized) {
		throw RuntimeError("DistortedPowerCorrelation::getCorrelationMultipole: not initialized.");
	}
	if(ell < 0 || ell > _ellMax || (_symmetric && (ell%2))) {
		throw RuntimeError("DistortedPowerCorrelation::getCorrelationMultipole: invalid ell.");
	}
	if(r < _rgrid.front() || r > _rgrid.back()) {
		throw RuntimeError("DistortedPowerCorrelation::getCorrelationMultipole: r out of range.");
	}
	int idx = _symmetric ? ell/2 : ell;
	return (*_interpolator[idx])(r);
}

double local::DistortedPowerCorrelation::getCorrelation(double r, double mu) const {
	if(!_initialized) {
		throw RuntimeError("DistortedPowerCorrelation::getCorrelation: not initialized.");
	}
	if(mu < 0 || mu > 1) {
		throw RuntimeError("DistortedPowerCorrelation::getPower: expected 0 <= mu <= 1.");
	}
	double result(0);
	int dell = _symmetric ? 2 : 1;
	for(int ell = 0; ell <= _ellMax; ell += dell) {
		result += getCorrelationMultipole(r,ell)*legendreP(ell,mu);
	}
	return result;
}

void local::DistortedPowerCorrelation::initialize() {
	// Loop over multipoles
	int dell = _symmetric ? 2 : 1;
	for(int ell = 0; ell <= _ellMax; ell += dell) {
		int idx(ell/dell);
		// Build a function object that evaluates this multipole for arbitrary k
		likely::GenericFunctionPtr fOfKPtr(
			new likely::GenericFunction(boost::bind(
				&DistortedPowerCorrelation::getPowerMultipole,this,_1,ell)));
		std::cout << "initializing ell = " << ell << std::endl;
		_transformer[idx]->initialize(fOfKPtr,_xiMoments[idx]);
	}
	_initialized = true;
}

bool local::DistortedPowerCorrelation::transform(bool bypassTerminationTest) const {
	// Loop over multipoles
	int dell = _symmetric ? 2 : 1;
	for(int ell = 0; ell <= _ellMax; ell += dell) {
		int idx(ell/dell);
		// Build a function object that evaluates this multipole for arbitrary k
		likely::GenericFunctionPtr fOfKPtr(
			new likely::GenericFunction(boost::bind(
				&DistortedPowerCorrelation::getPowerMultipole,this,_1,ell)));
		std::cout << "transforming ell = " << ell << std::endl;
		_transformer[idx]->transform(fOfKPtr,_xiMoments[idx]);
		// (re)create the interpolator for this moment
		_interpolator[idx].reset(new likely::Interpolator(_rgrid,_xiMoments[idx],"cspline"));
	}
	return true;
}
