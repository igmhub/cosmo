// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/MultipoleTransform.h"
#include "cosmo/RuntimeError.h"

#include <boost/math/special_functions/gamma.hpp>

#include <cmath>

namespace local = cosmo;

local::MultipoleTransform::MultipoleTransform(likely::GenericFunctionPtr func,
Type type, int ell, double vmin, double vmax, double veps, int minSamplesPerDecade)
{
	// Input parameter validation
	if(type != SphericalBessel && type != Hankel) {
		throw RuntimeError("MultipoleTransform: invalid type.");
	}
	if(ell < 0) {
		throw RuntimeError("MultipoleTransform: expected ell >= 0.");
	}
	if(vmin >= vmax) {
		throw RuntimeError("MultipoleTransform: expected vmin < vmax.");
	}
	if(veps <= 0) {
		throw RuntimeError("MultipoleTransform: expected veps > 0.");
	}
	double pi(atan2(0,-1));
	double alpha, uv0, s0;
	if(type == SphericalBessel) {
		// Calculate alpha and uv0 of eqn (1.6)
		alpha = 0.5*(1-ell);
		double gammaEll32 = boost::math::tgamma(ell+1.5);
		uv0 = 2*std::pow(gammaEll32/std::sqrt(pi),1./(ell+1));
		// Calculate s0 of eqn (1.8)
		s0 = 2./(ell+1);
	}
	else {
		// Calculate alpha and uv0 of eqn (2.4)
		alpha = 0.25*(1-2*ell);
		double gammaEll1 = boost::math::tgamma(ell+1);
		uv0 = 2*std::pow(gammaEll1/std::sqrt(pi),1./(ell+0.5));
		// Calculate s0 of eqn (2.6)
		s0 = 4./(2*ell+1);
	}
	// Calculate c of eqn (3.4)
	double c = pi/uv0;
	// Calculate eps from veps using the approx of eqn (3.7)
	double L0 = veps/c;
	if(L0 > 0.35) {
		throw RuntimeError("MultipoleTransform: veps to large for eqn (3.9) approx.");
	}
	double L1 = std::log(L0), L2 = std::log(-L1), L1sq(L1*L1);
	double arg = 6*L1sq*L1sq + 6*L1sq*L2*(L1+1) - 3*L1*L2*(L2-2) + L2*(2*L2*L2-9*L2+6);
	double eps = std::pow(-L0/(6*L1sq*L1)*arg,1./s0);
	// Calculate Y of eqn (3.3)
	double Y = uv0/(2*pi)*std::pow(eps,-s0);
	// Calculate delta of eqn (3.2)
	if(type == SphericalBessel) {
		arg = std::ceil(Y)/Y;
	}
	else {
		arg = (std::ceil(1./8.+Y/4.) - 1./8.)/Y;
	}
	double delta = std::log(arg);
	// Calculate sN of eqn (3.2)
	double sN = -s0*std::log(eps) + delta;
	// Calculate dsmax of eqn (3.4)
	double dsmax = c*std::pow(eps,s0);
	// Calculate the ds value corresponding to the min required number of samples per
	// decade, and use the smaller of these.
	double dsmax2 = std::log(10)/minSamplesPerDecade;
	if(dsmax2 < dsmax) dsmax = dsmax2;

	// Calculate the geometric mean of the target v range of eqn (3.9)
	double v0 = std::sqrt(vmin*vmax);
	// Calculate the corresponding u0
	double u0 = uv0/v0;

}

local::MultipoleTransform::~MultipoleTransform() { }
