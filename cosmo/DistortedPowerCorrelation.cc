// Created 20-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/DistortedPowerCorrelation.h"
#include "cosmo/AdaptiveMultipoleTransform.h"

#include "cosmo/RuntimeError.h"

namespace local = cosmo;

local::DistortedPowerCorrelation::DistortedPowerCorrelation(
TabulatedPowerCPtr power, RMuFunctionCPtr distortion, double rmin, double rmax, int nr,
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
	_moments.reserve(symmetric ? ellMax/2 : ellMax);
	for(int ell = 0; ell < ellMax; ell += dell) {
		AdaptiveMultipoleTransformPtr amt(new AdaptiveMultipoleTransform(
			MultipoleTransform::SphericalBessel,ell,_rgrid,relerr,abserr,abspow));
		_moments.push_back(amt);
	}
}

local::DistortedPowerCorrelation::~DistortedPowerCorrelation() { }

void local::DistortedPowerCorrelation::initialize() {

}
