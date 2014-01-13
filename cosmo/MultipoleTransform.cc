// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/MultipoleTransform.h"
#include "cosmo/RuntimeError.h"

#include "config.h"
#ifdef HAVE_LIBFFTW3F
#include "fftw3.h"
//#define FFTW(X) fftw_ ## X // double transforms
#define FFTW(X) fftwf_ ## X // float transforms
//#define FFTW(X) fftwl_ ## X // long double transforms
typedef float FftwReal;
#endif

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include <cmath>
#include <cstdlib> // for abs(int)
#include <vector>

namespace local = cosmo;

namespace cosmo {
    struct MultipoleTransform::Implementation {
#ifdef HAVE_LIBFFTW3F
        FFTW(complex) *fdata,*gdata;
        FFTW(plan) fplan,gplan,fgplan;
#endif
    };
} // cosmo::

local::MultipoleTransform::MultipoleTransform(Type type, int ell,
double vmin, double vmax, double veps, Strategy strategy,
int minSamplesPerCycle, int minSamplesPerDecade) :
_type(type),_minSamplesPerCycle(minSamplesPerCycle),
_pimpl(new Implementation())
{
#ifndef HAVE_LIBFFTW3F
	throw RuntimeError("MultipoleTransform: library not built with fftw3f support.");
#endif
	// Input parameter validation
	if(_type != SphericalBessel && _type != Hankel) {
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
	if(veps == 0) {
		throw RuntimeError("MultipoleTransform: expected veps != 0.");
	}
	if(minSamplesPerCycle <= 0) {
		throw RuntimeError("MultipoleTransform: expected minSamplesPerCycle > 0.");
	}
	if(minSamplesPerDecade <= 0) {
		throw RuntimeError("MultipoleTransform: expected minSamplesPerDecade > 0.");
	}
	double pi(atan2(0,-1));
	double alpha, uv0, s0;
	if(_type == SphericalBessel) {
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
	double arg, c = 2*pi/minSamplesPerCycle/uv0;
	if(veps > 0) {
		// Calculate eps from veps using the approx of eqn (3.7)
		double L0 = veps/c;
		if(L0 > 0.35) {
			throw RuntimeError("MultipoleTransform: veps to large for eqn (3.9) approx.");
		}
		double L1 = std::log(L0), L2 = std::log(-L1), L1sq(L1*L1);
		arg = 6*L1sq*L1sq + 6*L1sq*L2*(L1+1) - 3*L1*L2*(L2-2) + L2*(2*L2*L2-9*L2+6);
		_eps = std::pow(-L0/(6*L1sq*L1)*arg,1./s0);
	}
	else {
		_eps = -veps;
	}
	// Calculate Y of eqn (3.3)
	double Y = uv0/(2*pi)*std::pow(_eps,-s0);
	// Calculate delta of eqn (3.2)
	if(_type == SphericalBessel) {
		arg = std::ceil(Y)/Y;
	}
	else {
		arg = (std::ceil(1./8.+Y/4.) - 1./8.)/Y;
	}
	double delta = std::log(arg);
	// Calculate sN of eqn (3.2)
	double sN = -s0*std::log(_eps) + delta;
	// Calculate dsmax of eqn (3.4)
	double dsmax = c*std::pow(_eps,s0);
	// Calculate the ds value corresponding to the min required number of samples per
	// decade, and use the smaller of these.
	double dsmaxAlt = std::log(10)/minSamplesPerDecade;
	if(dsmaxAlt < dsmax) dsmax = dsmaxAlt;
	// Calculate Nf and ds of eqn (3.5)
	_Nf = (int)std::ceil(sN/dsmax);
	double ds = sN/_Nf;
	// Calculate the geometric mean of the target v range of eqn (3.9)
	double v0 = std::sqrt(vmin*vmax);
	// Calculate the corresponding u0 and its powers
	double u0 = uv0/v0;
	double u02 = u0*u0, u03 = u0*u02;
	// Calculate Ng of eqn (3.1)
	int Ng = (int)std::ceil(std::log(vmax/vmin)/(2*ds));
	// Tabulate f(s) of eqn (1.4) or (2.2)
	int Ntot = _Nf + Ng;
#ifdef HAVE_LIBFFTW3F
	// Allocate SIMD aligned arrays using FFTW's allocator
	_pimpl->fdata = (FFTW(complex)*)FFTW(malloc)(sizeof(FFTW(complex))*2*Ntot);
	_pimpl->gdata = (FFTW(complex)*)FFTW(malloc)(sizeof(FFTW(complex))*2*Ntot);
	// Build plans for doing transforms in place. FFTW is supposed to be smart
	// about re-using plans with the same dimensions but different arrays.
	int flags = (strategy == EstimatePlan) ? FFTW_ESTIMATE : FFTW_MEASURE;
	_pimpl->fplan = FFTW(plan_dft_1d)(2*Ntot,_pimpl->fdata,_pimpl->fdata,
		FFTW_FORWARD,flags);
	_pimpl->gplan = FFTW(plan_dft_1d)(2*Ntot,_pimpl->gdata,_pimpl->gdata,
		FFTW_FORWARD,flags);
	_pimpl->fgplan = FFTW(plan_dft_1d)(2*Ntot,_pimpl->gdata,_pimpl->gdata,
		FFTW_BACKWARD,flags);
	for(int m = 0; m < 2*Ntot; ++m) {
		int n = m;
		if(n >= Ntot) n -= 2*Ntot;
		if(std::abs(n) > _Nf) {
			_pimpl->fdata[m][0] = _pimpl->fdata[m][1] = 0.;
		}
		else {
			double bessel,s = n*ds;
			arg = uv0*std::exp(s);
			if(_type == SphericalBessel) {
				bessel = boost::math::sph_bessel(ell,arg);
			}
			else {
				bessel = boost::math::cyl_bessel_j(ell,arg);
			}
			_pimpl->fdata[m][0] = std::exp(alpha*s)*bessel*ds;
			_pimpl->fdata[m][1] = 0.;
		}
	}
	// Calculate the Fourier transform of fdata
	FFTW(execute)(_pimpl->fplan);
#endif
	// Tabulate the u values where func(u) should be evaluated, the
	// coefficients needed to rescale func(u(s)) to g(s), the v values
	// for the convolution result, and the scale factors for the result.
	_ugrid.reserve(2*Ntot);
	_vgrid.reserve(2*Ntot);
	_coef.reserve(2*Ntot);
	_scale.reserve(2*Ntot);
	for(int n = -Ntot; n < Ntot; ++n) {
		double v, s = n*ds;
		_ugrid.push_back(u0*std::exp(-s));
		_vgrid.push_back(v = v0*std::exp(+s));
		if(_type == SphericalBessel) {
			_coef.push_back(ds*std::exp((3-alpha)*(-s))*u03);
		}
		else {
			_coef.push_back(ds*std::exp((2-alpha)*(-s))*u02);
		}
		_scale.push_back(std::pow(v/v0,-alpha)/ds);
	}
}

local::MultipoleTransform::MultipoleTransform(
MultipoleTransform const &other, int subsampling) :
_type(other._type),_pimpl(new Implementation())
{
}

local::MultipoleTransform::~MultipoleTransform() {
#ifdef HAVE_LIBFFTW3F
    FFTW(destroy_plan)(_pimpl->fplan);
    FFTW(destroy_plan)(_pimpl->gplan);
    FFTW(destroy_plan)(_pimpl->fgplan);
    FFTW(free)(_pimpl->fdata);
    FFTW(free)(_pimpl->gdata);
#endif
}

void local::MultipoleTransform::transform(std::vector<double> const &funcTable,
std::vector<double> &result) const {
#ifndef HAVE_LIBFFTW3F
	throw RuntimeError("MultipoleTransform: library not built with fftw3f support.");
#else
	int n(_ugrid.size());
	// (re)initialize result vector to have correct size, if necessary
	if(result.size() != n) std::vector<double>(n,0).swap(result);
	for(int m = 0; m < n; ++m) {
		_pimpl->gdata[m][0] = (FftwReal)(_coef[m]*funcTable[m]);
		_pimpl->gdata[m][1] = 0.;
	}
	// Calculate the Fourier transform of gdata
	FFTW(execute)(_pimpl->gplan);
	// Multiply the transforms of fdata and gdata, saving the result in gdata
	double norm(n);
	for(int m = 0; m < n; ++m) {
		double re1 = _pimpl->fdata[m][0], im1 = _pimpl->fdata[m][1];
		double re2 = _pimpl->gdata[m][0], im2 = _pimpl->gdata[m][1];
		_pimpl->gdata[m][0] = (FftwReal)((re1*re2 - im1*im2)/norm);
		_pimpl->gdata[m][1] = (FftwReal)((re1*im2 + re2*im1)/norm);
	}
	// Calculate the inverse Fourier transform that gives the convolution of
	// the original fdata and gdata, tabulated on vgrid.
	FFTW(execute)(_pimpl->fgplan);
	// Rescale and copy the results back to the vector provided.
	for(int m = 0; m < n; ++m) {
		result[m] = _scale[m]*_pimpl->gdata[m][0];
	}
#endif
}
