// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/MultipoleTransform.h"
#include "cosmo/RuntimeError.h"

#include <boost/math/special_functions/gamma.hpp>

#include <cmath>

namespace local = cosmo;

local::MultipoleTransform::MultipoleTransform(likely::GenericFunctionPtr func,
Type type, int ell, double vmin, double vmax, double eps)
{
	// Input parameter validation
	if(ell < 0) {
		throw RuntimeError("MultipoleTransform: expected ell >= 0.");
	}
	if(vmin >= vmax) {
		throw RuntimeError("MultipoleTransform: expected vmin < vmax.");
	}
	if(eps <= 0) {
		throw RuntimeError("MultipoleTransform: expected eps > 0.");
	}
	// Calculate alpha and uv0 of eqn (3.1)
	double alpha = 0.5*(1-ell);
	double pi(atan2(0,-1));
	double uv0 = 2*std::pow(boost::math::tgamma(ell+1.5)/std::sqrt(pi),1./(ell+1));
	// Calculate the geometric mean of the target v range of eqn (3.11)
	double v0 = std::sqrt(vmin*vmax);
	// Calculate the corresponding u0
	double u0 = uv0/v0;
}

local::MultipoleTransform::~MultipoleTransform() { }
