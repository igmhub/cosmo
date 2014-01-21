// Created 20-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/DistortedPowerCorrelation.h"

#include "cosmo/RuntimeError.h"

namespace local = cosmo;

local::DistortedPowerCorrelation::DistortedPowerCorrelation(
TabulatedPowerCPtr power, RMuFunctionCPtr distortion, int ellMax, double rmin, double rmax)
: _power(power), _distortion(distortion), _ellMax(ellMax), _rmin(rmin), _rmax(rmax)
{
	if(ellMax < 0) {
		throw RuntimeError("DistortedPowerCorrelation: expected ellMax >= 0.");
	}
	if(rmax <= rmin) {
		throw RuntimeError("DistortedPowerCorrelation: expected rmin < rmax.");
	}
	if(rmin <= 0) {
		throw RuntimeError("DistortedPowerCorrelation: expected rmin > 0.");
	}
}

local::DistortedPowerCorrelation::~DistortedPowerCorrelation() { }
