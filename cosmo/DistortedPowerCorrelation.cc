// Created 20-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/DistortedPowerCorrelation.h"
#include "cosmo/AdaptiveMultipoleTransform.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "cosmo/RuntimeError.h"

#include "likely/Interpolator.h"

#include "boost/foreach.hpp"
#include "boost/bind.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>

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
		// Use the same relerr for each ell and share abserr equally. These values will
		// be adjusted when initialize is called later.
		AdaptiveMultipoleTransformPtr amt(new AdaptiveMultipoleTransform(
			MultipoleTransform::SphericalBessel,ell,coef,_rgrid,relerr/10.,abserr/(2*nell),abspow));
		_transformer.push_back(amt);
		_xiMoments.push_back(std::vector<double>(nr,0.));
		_interpolator.push_back(likely::InterpolatorPtr());
	}
	// initialize vectors used to find biggest relative contributions
	std::vector<double>(nell).swap(_rbig);
	std::vector<double>(nell).swap(_mubig);
	std::vector<double>(nell).swap(_relbig);
}

local::DistortedPowerCorrelation::~DistortedPowerCorrelation() { }

double local::DistortedPowerCorrelation::getPower(double k, double mu) const {
	if(mu < -1 || mu > 1) {
		throw RuntimeError("DistortedPowerCorrelation::getPower: expected -1 <= mu <= 1.");
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

double local::DistortedPowerCorrelation::getSavedPowerMultipole(double k, int ell) const {
	if(ell < 0 || ell > _ellMax || (_symmetric && (ell%2))) {
		throw RuntimeError("DistortedPowerCorrelation::getPowerMultipole: invalid ell.");
	}
	return 0;
}

double local::DistortedPowerCorrelation::getCorrelationMultipole(double r, int ell) const {
	if(!isInitialized()) {
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
	if(!isInitialized()) {
		throw RuntimeError("DistortedPowerCorrelation::getCorrelation: not initialized.");
	}
	if(mu < -1 || mu > 1) {
		throw RuntimeError("DistortedPowerCorrelation::getPower: expected -1 <= mu <= 1.");
	}
	double result(0);
	int dell = _symmetric ? 2 : 1;
	for(int ell = 0; ell <= _ellMax; ell += dell) {
		result += getCorrelationMultipole(r,ell)*legendreP(ell,mu);
	}
	return result;
}

void local::DistortedPowerCorrelation::initialize(int nmu, int minSamplesPerDecade,
double margin, double vepsMax, double vepsMin, bool optimize) {
	if(nmu < 2) {
		throw RuntimeError("DistortedPowerCorrelation::initialize: expected nmu >= 2.");
	}
	if(minSamplesPerDecade < 0) {
		throw RuntimeError("DistortedPowerCorrelation::initialize: expected minSamplesPerDecade >= 0.");
	}
	if(margin < 1) {
		throw RuntimeError("DistortedPowerCorrelation::initialize: expected margin >= 1.");
	}
	if(vepsMax <= vepsMin) {
		throw RuntimeError("DistortedPowerCorrelation::initialize: expected vepsMax > vepsMin.");
	}
	if(vepsMin <= 0) {
		throw RuntimeError("DistortedPowerCorrelation::initialize: expected vepsMin > 0.");
	}
	// Loop over multipoles
	int dell = _symmetric ? 2 : 1;
	for(int ell = 0; ell <= _ellMax; ell += dell) {
		int idx(ell/dell);
		// Build a function object that evaluates this multipole for arbitrary k
		likely::GenericFunctionPtr fOfKPtr(
			new likely::GenericFunction(boost::bind(
				&DistortedPowerCorrelation::getPowerMultipole,this,_1,ell)));
		// Do not optimize now
		bool noOptimize(false);
		_transformer[idx]->initialize(fOfKPtr,_xiMoments[idx],minSamplesPerDecade,margin,
			vepsMax,vepsMin,noOptimize);
		// (re)create the interpolator for this moment
		_interpolator[idx].reset(new likely::Interpolator(_rgrid,_xiMoments[idx],"cspline"));
	}
	// Loop over our (r,mu) evaluation grid.
	double dmu = 2./dell/(nmu-1.);
	int nell = 1 + _ellMax/dell;
	std::vector<double> contribution(nell);
	// Clear our relative errors
	std::fill(_relbig.begin(),_relbig.end(),0.);
	BOOST_FOREACH(double r, _rgrid) {
		for(int i = 0; i < nmu; ++i) {
			double mu = 1. - i*dmu;
			// Loop over multipoles to calculate their relative contributions at (r,mu)
			double xisum;
			for(int ell = 0; ell <= _ellMax; ell += dell) {
				int idx = _symmetric ? ell/2 : ell;
				double term = (*_interpolator[idx])(r)*legendreP(ell,mu);
				contribution[idx] = term;
				xisum += term;
			}
			// Skip (r,mu) points where xi is essentially zero
			if(std::fabs(xisum) < _abserr*std::pow(r,_abspow)) continue;
			// Update the biggest relative contribution found so far for each multipole
			for(int idx = 0; idx < nell; ++idx) {
				double relfrac = std::fabs(contribution[idx]/xisum);
				if(relfrac > _relbig[idx]) {
					_rbig[idx] = r;
					_mubig[idx] = mu;
					_relbig[idx] = relfrac;
				}
			}
		}
	}
	// Reset our transformers using updated relerr specs
	for(int ell = 0; ell <= _ellMax; ell += dell) {
		int idx = _symmetric ? ell/2 : ell;
		double relerr = _relerr/nell/_relbig[idx];
		double abserr = _abserr/nell;
		double coef = multipoleTransformNormalization(ell,3,+1);
		AdaptiveMultipoleTransformPtr amt(new AdaptiveMultipoleTransform(
			MultipoleTransform::SphericalBessel,ell,coef,_rgrid,relerr,abserr,_abspow));
		_transformer[idx] = amt;
		// Build a function object that evaluates this multipole for arbitrary k
		likely::GenericFunctionPtr fOfKPtr(
			new likely::GenericFunction(boost::bind(
				&DistortedPowerCorrelation::getPowerMultipole,this,_1,ell)));
		// Initialize our new transformer (with optimization, if requested)
		_transformer[idx]->initialize(fOfKPtr,_xiMoments[idx],minSamplesPerDecade,margin,
			vepsMax,vepsMin,optimize);
		// (re)create the interpolator for this moment
		_interpolator[idx].reset(new likely::Interpolator(_rgrid,_xiMoments[idx],"cspline"));
	}
	_initialized = true;
}

bool local::DistortedPowerCorrelation::transform(bool bypassTerminationTest) const {
	bool accurate(true);
	// Loop over multipoles
	int dell = _symmetric ? 2 : 1;
	for(int ell = 0; ell <= _ellMax; ell += dell) {
		int idx(ell/dell);
		// Build a function object that evaluates this multipole for arbitrary k
		likely::GenericFunctionPtr fOfKPtr(
			new likely::GenericFunction(boost::bind(
				&DistortedPowerCorrelation::getPowerMultipole,this,_1,ell)));
		accurate &= _transformer[idx]->transform(fOfKPtr,_xiMoments[idx],bypassTerminationTest);
		// (re)create the interpolator for this moment
		_interpolator[idx].reset(new likely::Interpolator(_rgrid,_xiMoments[idx],"cspline"));
	}
	return accurate;
}

local::AdaptiveMultipoleTransformCPtr local::DistortedPowerCorrelation::getTransform(int ell) const {
	if(ell < 0 || ell > _ellMax || (_symmetric && (ell%2))) {
		throw RuntimeError("DistortedPowerCorrelation::getTransform: invalid ell.");
	}
	int dell = _symmetric ? 2 : 1;
	int idx(ell/dell);
	return _transformer[idx];
}

void local::DistortedPowerCorrelation::getBiggestContribution(int ell,
double &rbig, double &mubig, double &relbig) const {
	if(!isInitialized()) {
		throw RuntimeError("DistortedPowerCorrelation::getBiggestContribution: not initialized.");
	}
	if(ell < 0 || ell > _ellMax || (_symmetric && (ell%2))) {
		throw RuntimeError("DistortedPowerCorrelation::getBiggestContribution: invalid ell.");
	}
	int dell = _symmetric ? 2 : 1;
	int idx(ell/dell);
	rbig = _rbig[idx];
	mubig = _mubig[idx];
	relbig = _relbig[idx];
}

void local::DistortedPowerCorrelation::printToStream(std::ostream &out) const {
    double r,mu,rel;
    int dell = _symmetric ? 2 : 1;
    out << "xi(r,mu) interpolated at " << _rgrid.size() << " points covering r = ["
    	<< _rgrid.front() << ',' << _rgrid.back() << "] Mpc/h with even ell <= "
		<< _ellMax << std::endl;
    for(int ell = 0; ell <= _ellMax; ell += dell) {
        getBiggestContribution(ell,r,mu,rel);
        cosmo::AdaptiveMultipoleTransformCPtr amt = getTransform(ell);
        out << "initialized ell = " << ell << " adaptive transform:" << std::endl;
        out << "  relerr = " << amt->getRelErr() << " @(r=" << r << " Mpc/h,mu=" << mu
        	<< ",rel=" << rel << "), abserr = " << amt->getAbsErr() << " (abspow = "
            << amt->getAbsPow() << ")," << std::endl;
		out << "  veps = " << amt->getVEps() << ", kmin = " << amt->getUMin() << " h/Mpc, kmax = "
			<< amt->getUMax() << " h/Mpc, nk = " << amt->getNU() << " ("
			<< (int)std::floor(amt->getUSamplesPerDecade()) << " samples/decade)" << std::endl;
    }	
}
