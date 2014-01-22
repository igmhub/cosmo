// Created 20-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/DistortedPowerCorrelation.h"
#include "cosmo/AdaptiveMultipoleTransform.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "cosmo/RuntimeError.h"

#include "boost/foreach.hpp"
#include "boost/bind.hpp"

namespace local = cosmo;

local::DistortedPowerCorrelation::DistortedPowerCorrelation(likely::GenericFunctionPtr power,
RMuFunctionCPtr distortion, double rmin, double rmax, int nr,
int ellMax, bool symmetric, double relerr, double abserr, double abspow)
: _power(power), _distortion(distortion), _ellMax(ellMax), _symmetric(symmetric),
_relerr(relerr), _abserr(abserr), _abspow(abspow)
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
	for(int ell = 0; ell <= ellMax; ell += dell) {
		// use the same relerr for each ell and share abserr equally
		AdaptiveMultipoleTransformPtr amt(new AdaptiveMultipoleTransform(
			MultipoleTransform::SphericalBessel,ell,_rgrid,relerr,abserr/nell,abspow));
		_transformer.push_back(amt);
		_xiMoments.push_back(std::vector<double>(nr,0.));
	}
}

local::DistortedPowerCorrelation::~DistortedPowerCorrelation() { }

double local::DistortedPowerCorrelation::getPower(double k, double mu) const {
	return (*_power)(k)*(*_distortion)(k,mu);
}

double local::DistortedPowerCorrelation::getPowerMultipole(double k, int ell) const {
	likely::GenericFunctionPtr fOfMuPtr(
		new likely::GenericFunction(boost::bind(
			&DistortedPowerCorrelation::getPower,this,k,_1)));
	return getMultipole(fOfMuPtr, ell);
}

void local::DistortedPowerCorrelation::initialize() {
	// Loop over multipoles
	int dell = _symmetric ? 2 : 1;
	for(int ell = 0; ell <= _ellMax; ell += dell) {
		// Build a function object that evaluates this multipole for arbitrary k
		likely::GenericFunctionPtr fOfKPtr(
			new likely::GenericFunction(boost::bind(
				&DistortedPowerCorrelation::getPowerMultipole,this,_1,ell)));
		std::cout << "initializing ell = " << ell << std::endl;
		_transformer[ell/dell]->initialize(fOfKPtr,_xiMoments[ell/dell]);
	}
}

bool local::DistortedPowerCorrelation::transform(bool bypassTerminationTest) const {
	// Loop over multipoles
	int dell = _symmetric ? 2 : 1;
	for(int ell = 0; ell <= _ellMax; ell += dell) {
		// Build a function object that evaluates this multipole for arbitrary k
		likely::GenericFunctionPtr fOfKPtr(
			new likely::GenericFunction(boost::bind(
				&DistortedPowerCorrelation::getPowerMultipole,this,_1,ell)));
		std::cout << "transforming ell = " << ell << std::endl;
		_transformer[ell/dell]->transform(fOfKPtr,_xiMoments[ell/dell]);
	}
	return true;
}
