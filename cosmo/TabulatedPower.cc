// Created 13-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/TabulatedPower.h"
#include "cosmo/RuntimeError.h"

#include "likely/Interpolator.h"

#include <cmath>
#include <iostream>
#include <fstream>

namespace local = cosmo;

namespace cosmo {
	struct TabulatedPower::PowerLawExtrapolator {
		// Performs a power-law extrapolation P(k) = c*k^p using parameters
		// c,p determined from two points (k1,P1=P(k1)) and (k2,P2=P(k2)).
		PowerLawExtrapolator(double k1,double P1,double k2,double P2) :
		a(std::log(P2/P1)/std::log(k2/k1)), c(P1/std::pow(k1,a)) { }
		double operator()(double k) const { return c*std::pow(k,a); }
		double a,c;
	};
} // cosmo::

local::TabulatedPower::TabulatedPower(
std::vector<double> const &k, std::vector<double> const &Pk,
bool extrapolateBelow, bool extrapolateAbove, double maxRelError, bool verbose) :
_maxRelError(maxRelError)
{
	// Check that input vectors have the same size
	if(k.size() != Pk.size()) {
		throw RuntimeError("TabulatedPower: input vectors have different sizes.");
	}
	if(k.size() < 3 && (extrapolateBelow || extrapolateAbove)) {
		throw RuntimeError("TabulatedPower: need at least 3 points for extrapolation.");
	}
	// Convert k to log(k)
	std::vector<double> logk;
	logk.reserve(k.size());
	double lastk(0);
	for(std::vector<double>::const_iterator iter = k.begin(); iter != k.end(); ++iter) {
		double k(*iter);
		if(k <= lastk) {
			throw RuntimeError("TabulatedPower: invalid input k vector.");
		}
		logk.push_back(std::log(k));
	}
	// Remember our interpolation limits
	_kmin = k.front();
	_kmax = k.back();
	if(verbose) {
		std::cout << "TabulatedPower: using " << k.size() << " points covering "
			<< _kmin << " <= k <= " << _kmax << std::endl;
	}
	// Build a spline interpolator in log(k) and P(k)
	_interpolator.reset(new likely::Interpolator(logk,Pk,"cspline"));
	// Estimate a power law for extrapolating below kmin, if requested
	if(extrapolateBelow) {
		_extrapolateBelow.reset(new PowerLawExtrapolator(k[0],Pk[0],k[2],Pk[2]));
		// Check how well the extrapolation does at k[1]
		double P1 = (*_extrapolateBelow)(k[1]);
		double relerr = std::fabs(P1/Pk[1]-1.);
		if(verbose) {
			std::cout << "TabulatedPower: relative error for extrapolation below is "
				<< relerr << std::endl;
		}
		if(relerr > maxRelError) {
			throw RuntimeError("TabulatedPower: cannot reliably extrapolate below kmin.");
		}
	}
	// Estimate a power law for extrapolating above kmax, if requested
	if(extrapolateAbove) {
		int n = k.size();
		_extrapolateAbove.reset(new PowerLawExtrapolator(k[n-3],Pk[n-3],k[n-1],Pk[n-1]));
		// Check how well the extrapolation does at k[n-2]
		double Pn2 = (*_extrapolateAbove)(k[n-2]);
		double relerr = std::fabs(Pn2/Pk[n-2]-1.);
		if(verbose) {
			std::cout << "TabulatedPower: relative error for extrapolation above is "
				<< relerr << std::endl;
		}
		if(relerr > maxRelError) {
			throw RuntimeError("TabulatedPower: cannot reliably extrapolate above kmax.");
		}
	}
}

local::TabulatedPower::~TabulatedPower() { }

double local::TabulatedPower::operator()(double k) const {
	// Return 0 without complaining if k <= 0
	if(k <= 0) return 0;
	// Extrapolate if enabled.
	if(k < _kmin) {
		if(_extrapolateBelow) {
			return (*_extrapolateBelow)(k);
		}
		else {
			throw RuntimeError("TabulatedPower: extrapolation below kmin not enabled.");
		}
	}
	else if(k > _kmax) {
		if(_extrapolateAbove) {
			return (*_extrapolateAbove)(k);
		}
		else {
			throw RuntimeError("TabulatedPower: extrapolation above kmax not enabled.");
		}
	}
	else {
		// We must have kmin <= k <= kmax so interpolate in log(k)
		double logk = std::log(k);
		return (*_interpolator)(logk);
	}
}

local::TabulatedPowerCPtr local::TabulatedPower::createDelta(TabulatedPowerCPtr other) const {
	likely::Interpolator::CoordinateValues logkGrid = _interpolator->getXGrid();
	likely::Interpolator::CoordinateValues deltaGrid = _interpolator->getYGrid();
	likely::Interpolator::CoordinateValues kGrid;
	int n(logkGrid.size());
	kGrid.reserve(n);
	for(int i = 0; i < n; ++i) {
		double k = std::exp(logkGrid[i]);
		kGrid.push_back(k);
		deltaGrid[i] -= (*other)(k);
	}
	bool extrapolateBelow(!!_extrapolateBelow), extrapolateAbove(!!_extrapolateAbove);
	TabulatedPowerCPtr delta(new TabulatedPower(kGrid,deltaGrid,
		extrapolateBelow,extrapolateAbove,_maxRelError));
	return delta;
}

local::TabulatedPowerCPtr local::createTabulatedPower(std::string const &filename,
bool extrapolateBelow, bool extrapolateAbove, double maxRelError, bool verbose)
{
	std::vector<std::vector<double> > columns(2);
    std::ifstream input(filename.c_str());
    likely::readVectors(input, columns);
    TabulatedPowerCPtr power(new TabulatedPower(columns[0],columns[1],
    	extrapolateBelow,extrapolateAbove,maxRelError,verbose));
    return power;
}
