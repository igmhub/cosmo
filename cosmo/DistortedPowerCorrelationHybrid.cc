// Created 11-May-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#include "cosmo/DistortedPowerCorrelationHybrid.h"
#include "cosmo/RuntimeError.h"

#include "likely/BiCubicInterpolator.h"
#include "likely/Integrator.h"
#include "likely/Interpolator.h"

#include "boost/bind.hpp"
#include <boost/math/special_functions/bessel.hpp>

#include <cmath>
#include <algorithm>
#include <iostream>

#include "config.h"
#ifdef HAVE_LIBFFTW3F
#include "fftw3.h"
#define FFTW(X) fftwf_ ## X // prefix identifier (float transform)
#endif

namespace local = cosmo;

namespace cosmo {
    struct DistortedPowerCorrelationHybrid::Implementation {
#ifdef HAVE_LIBFFTW3F
        FFTW(complex) *data;
        FFTW(plan) plan;
#endif
    };
}

local::DistortedPowerCorrelationHybrid::DistortedPowerCorrelationHybrid(likely::GenericFunctionPtr power,
KMuPkFunctionCPtr distortion, double kxmin, double kxmax, int nx, double spacing, int ny, double rmax,
double epsAbs)
: _power(power), _distortion(distortion), _kxmin(0), _kxmax(kxmax), _nx(nx), _spacing(spacing), _ny(ny),
_rmax(rmax), _epsAbs(epsAbs), _pimpl(new Implementation())
{	
	// Input parameter validation.
	if(kxmax < kxmin) {
		throw RuntimeError("DistortedPowerCorrelationHybrid: expected kxmax > kxmin");
	}
	if(kxmin <= 0) {
		throw RuntimeError("DistortedPowerCorrelationHybrid: expected kxmin > 0.");
	}
	if(spacing <= 0) {
		throw RuntimeError("DistortedPowerCorrelationHybrid: invalid grid spacing.");
	}
	if(nx <= 0 || ny <= 0) {
		throw RuntimeError("DistortedPowerCorrelationHybrid: invalid grid size.");
	}
	if(rmax <= 0) {
		throw RuntimeError("DistortedPowerCorrelationHybrid: expected rmax > 0.");
	}
	if(rmax > ny*spacing/2.) {
		throw RuntimeError("DistortedPowerCorrelationHybrid: expected rmax < ny*spacing/2");
	}
#ifdef HAVE_LIBFFTW3F
	// Allocate data array for 1D Fourier transform.
	_pimpl->data = (FFTW(complex)*) FFTW(malloc)(sizeof(FFTW(complex)) * ny);
#else
    throw RuntimeError("DistortedPowerCorrelationHybrid: package not built with FFTW3.");
#endif
	_twopi = 8*std::atan(1);
    // Initialize the k-space grid that will be used for evaluating the power spectrum.
    // Note that the grid is centered on (kxmin,0) and grid points for the y axis wrap around
    // in a specific order.
    //double dkx = std::pow(kxmax/kxmin,1./(nx-1.));
    _dkx = kxmax/(nx-1);
    double dky = _twopi/(ny*spacing);
	_kxgrid.reserve(nx);
	_kygrid.reserve(ny);
	for(int ix = 0; ix < nx; ++ix){
        //double kx = kxmin*std::pow(dkx,ix);
        double kx = ix*_dkx;
        _kxgrid.push_back(kx);
    }
    for(int iy = 0; iy < ny; ++iy){
        double ky = (iy > ny/2 ? iy-ny : iy)*dky;
        _kygrid.push_back(ky);
    }
    // Calculate normalization.
    _norm = ny*spacing;
    // Initialize the array that will be used to store the results of the 1D FFTs.
    _ktf.reset(new double[nx*(ny/2+1)]);
    // Initialize the r-space grid that will be used for evaluating the correlation function.
    _nr = (int)std::ceil(rmax/spacing)+1;
    _rgrid.reserve(_nr);
	for(int i = 0; i < _nr; ++i){
        _rgrid.push_back(i*spacing);
    }
    // Initialize the array that will be used for bicubic interpolation of the correlation function.
    _xi.reset(new double[_nr*_nr]);
}

local::DistortedPowerCorrelationHybrid::~DistortedPowerCorrelationHybrid() {
#ifdef HAVE_LIBFFTW3F
    if(0 != _pimpl->data) {
        FFTW(destroy_plan)(_pimpl->plan);
        FFTW(free)(_pimpl->data);
    }
#endif
}

double local::DistortedPowerCorrelationHybrid::getPower(double k, double mu) const {
	if(mu < -1 || mu > 1) {
		throw RuntimeError("DistortedPowerCorrelationHybrid::getPower: expected -1 <= mu <= 1.");
	}
	if(k < 0 ) {
		throw RuntimeError("DistortedPowerCorrelationHybrid::getPower: expected k >= 0.");
	}
	return (*_power)(k)*(*_distortion)(k,mu,(*_power)(k));
}

double local::DistortedPowerCorrelationHybrid::getCorrelation(double r, double mu) const {
	if(mu < -1 || mu > 1) {
		throw RuntimeError("DistortedPowerCorrelationHybrid::getCorrelation: expected -1 <= mu <= 1.");
	}
	double rpar = r*std::fabs(mu);
	double rperp = r*std::sqrt(1-mu*mu);
	if(r < 0 || rperp > _rmax || rpar > _rmax) {
		throw RuntimeError("DistortedPowerCorrelationHybrid::getCorrelation: r out of range.");
	}
	return (*_xiInterpolator)(rperp,rpar);
}

double local::DistortedPowerCorrelationHybrid::getKTransform(double ry, double kx) const {
	if(kx < _kxmin || kx > _kxmax) {
		throw RuntimeError("DistortedPowerCorrelationHybrid::getKTransform: expected kxmin <= kx <= kxmax");
	}
	if(ry < 0 || ry > _rmax) {
		throw RuntimeError("DistortedPowerCorrelationHybrid::getKTransform: expected 0 <= ry <= rmax");
	}
	return (*_ktransformInterpolator)(ry,kx);
}

void local::DistortedPowerCorrelationHybrid::transform() {
    // Perform a series of 1D Fourier transforms
    ktransform();
    // Perform a series of 1D integrals.
    likely::Integrator::IntegrandPtr integrand(new likely::Integrator::Integrand(
        boost::bind(&DistortedPowerCorrelationHybrid::_transverseIntegrand,this,_1)));
    likely::Integrator integrator(integrand,_epsAbs,0);
    likely::Interpolator::CoordinateValues ktfrow(_nx);
    for(int iy = 0; iy < _nr; ++iy){
        for(int i = 0; i < _nx; ++i){
            std::size_t ind2(iy+(_ny/2+1)*i);
            ktfrow[i] = _ktf[ind2];
        }
        _ktfInterpolator.reset(new likely::Interpolator(_kxgrid,ktfrow,"cspline"));
        for(int ix = 0; ix < _nr; ++ix){
            _rx = _rgrid[ix];
            std::size_t ind(ix+_nr*iy);
            _xi[ind] = integrator.integrateSmooth(_kxmin,_kxmax)/_twopi;
            std::cout << _rgrid[iy] << " " << _rx << " " << _xi[ind] << std::endl;
        }
    }
    // Create the bicubic interpolator.
	_xiInterpolator = new likely::BiCubicInterpolator(likely::BiCubicInterpolator::DataPlane(_xi),_spacing,_nr);
}

double local::DistortedPowerCorrelationHybrid::_transverseIntegrand(double kx) const {
    return kx*boost::math::cyl_bessel_j(0,kx*_rx)*(*_ktfInterpolator)(kx);
}

void local::DistortedPowerCorrelationHybrid::ktransform() {
#ifdef HAVE_LIBFFTW3F
	for(int ix = 0; ix < _nx; ++ix){
	    // Clean up any previous plan.
        if(_pimpl->data) FFTW(destroy_plan)(_pimpl->plan);
        // Create a new plan for an in-place transform.
        _pimpl->plan = FFTW(plan_dft_1d)(_ny,_pimpl->data,_pimpl->data,FFTW_BACKWARD,FFTW_ESTIMATE);
	    // Evaluate the power spectrum at each grid point.
	    for(int iy = 0; iy < _ny; ++iy){
		    std::size_t index(iy);
			double ksq = _kxgrid[ix]*_kxgrid[ix] + _kygrid[iy]*_kygrid[iy];
            double k = std::sqrt(ksq);
            if(k==0) {
            	_pimpl->data[index][0] = 0;
            	_pimpl->data[index][1] = 0;
            }
            if(k>0) {
            	double mu = _kygrid[iy]/k;
            	_pimpl->data[index][0] = getPower(k,mu);
            	_pimpl->data[index][1] = 0;
            }
        }
        // Execute 1D FFT to r-space.
	    FFTW(execute)(_pimpl->plan);
	    // Extract the transform result for the positive quadrant.
	    for(int iy = 0; iy < _ny/2+1; ++iy) {
	        std::size_t ind(iy);
	        std::size_t ind2(iy+(_ny/2+1)*ix);
	        _ktf[ind2] = (double)_pimpl->data[ind][0]/_norm;
	    }
    }
    // Create the bicubic interpolator.
	_ktransformInterpolator = new likely::BiCubicInterpolator(likely::BiCubicInterpolator::DataPlane(_ktf),_spacing,_ny/2+1,_nx,_dkx);
#endif
}

std::size_t local::DistortedPowerCorrelationHybrid::getMemorySize() const {
    return sizeof(*this) + (std::size_t)_nx*_ny*8 + (std::size_t)_nr*_nr*8;
}
