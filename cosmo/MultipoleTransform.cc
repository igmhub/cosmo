// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/MultipoleTransform.h"
#include "cosmo/RuntimeError.h"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include <cmath>
#include <cstdlib> // for abs(int)
#include <vector>

namespace local = cosmo;

local::MultipoleTransform::MultipoleTransform(Type type, int ell,
double vmin, double vmax, double veps, int minSamplesPerDecade)
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
	if(vmin <= 0) {
		throw RuntimeError("MultipoleTransform: expected vmin > 0.");
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
	double dsmaxAlt = std::log(10)/minSamplesPerDecade;
	if(dsmaxAlt < dsmax) dsmax = dsmaxAlt;
	// Calculate Nf and ds of eqn (3.5)
	int Nf = (int)std::ceil(sN/dsmax);
	double ds = sN/Nf;
	// Calculate the geometric mean of the target v range of eqn (3.9)
	double v0 = std::sqrt(vmin*vmax);
	// Calculate the corresponding u0 and its powers
	double u0 = uv0/v0;
	double u02 = u0*u0, u03 = u0*u02;
	// Calculate Ng of eqn (3.1)
	int Ng = (int)std::ceil(std::log(vmax/vmin)/(2*ds));
	// Tabulate f(s) of eqn (1.4) or (2.2)
	int Ntot = Nf + Ng;
	std::vector<double> fdata(2*Ntot,0.);
	for(int m = 0; m < Ntot; ++m) {
		int n = m;
		if(n >= Ntot) n -= 2*Ntot;
		if(std::abs(n) > Nf) continue;
		double bessel,s = n*ds;
		arg = uv0*std::exp(s);
		if(type == SphericalBessel) {
			bessel = boost::math::sph_bessel(ell,arg);
		}
		else {
			bessel = boost::math::cyl_bessel_j(ell,arg);
		}
		fdata[m] = std::exp(alpha*s)*bessel*ds;
	}
	// Calculate the Fourier transform of fdata
	// ...
	// Tabulate the u values where func(u) should be evaluated, the
	// coefficients needed to rescale func(u(s)) to g(s), the v values
	// for the convolution result, and the scale factors for the result.
	std::vector<double> ugrid(2*Ntot), vgrid(2*Ntot), coef(2*Ntot), scale(2*Ntot);
	for(int n = -Ntot; n < Ntot; ++n) {
		double s = n*ds;
		ugrid[n] = u0*std::exp(-s);
		vgrid[n] = v0*std::exp(+s);
		if(type == SphericalBessel) {
			coef[n] = std::exp((3-alpha)*(-s))*u03;
		}
		else {
			coef[n] = std::exp((2-alpha)*(-s))*u02;
		}
		scale[n] = std::pow(vgrid[n]/v0,-alpha)/ds;
	}
}

local::MultipoleTransform::~MultipoleTransform() { }

void local::MultipoleTransform::transform(likely::GenericFunctionPtr func,
std::vector<double> const &vgrid, std::vector<double> &result) const {
}
