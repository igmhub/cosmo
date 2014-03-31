// Created 20-Mar-2014 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#include "cosmo/DistortedPowerCorrelationFft.h"
#include "cosmo/RuntimeError.h"

#include "likely/BiCubicInterpolator.h"

#include <cmath>
#include <algorithm>
#include <iostream>

#include "config.h"
#ifdef HAVE_LIBFFTW3
#include "fftw3.h"
#define FFTW(X) fftw_ ## X // prefix identifier (double precision transform)
#endif

namespace local = cosmo;

namespace cosmo {
    struct DistortedPowerCorrelationFft::Implementation {
#ifdef HAVE_LIBFFTW3
        FFTW(complex) *data;
        FFTW(plan) plan;
#endif
    };
}

local::DistortedPowerCorrelationFft::DistortedPowerCorrelationFft(likely::GenericFunctionPtr power,
RMuFunctionCPtr distortion, double spacing, int nx, int ny, int nz)
: _power(power), _distortion(distortion), _spacing(spacing), _nx(nx), _ny(ny), _nz(nz), _pimpl(new Implementation())
{	
#ifdef HAVE_LIBFFTW3
	// allocate data array
	_pimpl->data = (FFTW(complex)*) FFTW(malloc)(sizeof(FFTW(complex)) * nx*ny*nz);
#else
    throw RuntimeError("DistortedPowerCorrelationFft: package not built with FFTW3.");
#endif
	double twopi(8*std::atan(1));
    // initialize the k space grid that will be used for evaluating the power spectrum
    double dkx(twopi/(nx*spacing)), dky(twopi/(ny*spacing)), dkz(twopi/(nz*spacing));
    double kmin = 0.000001;
	_kxgrid.reserve(nx);
	_kygrid.reserve(ny);
	_kzgrid.reserve(nz);
	for(int ix = 0; ix < nx; ++ix){
        double kx = ix*dkx + kmin;
        _kxgrid.push_back(kx);
    }
    for(int iy = 0; iy < ny; ++iy){
        double ky = iy*dky + kmin;
        _kygrid.push_back(ky);
    }
    for(int iz = 0; iz < nz; ++iz){
        double kz = iz*dkz + kmin;
        _kzgrid.push_back(kz);
    }
    // calculate normalisation
    _norm = nx*ny*nz*spacing*spacing*spacing;
    // initialize the array that will be used for bicubic interpolation
    _xi.reset(new double[nx*ny]);
}

local::DistortedPowerCorrelationFft::~DistortedPowerCorrelationFft() {
#ifdef HAVE_LIBFFTW3
    if(0 != _pimpl->data) {
        FFTW(destroy_plan)(_pimpl->plan);
        FFTW(free)(_pimpl->data);
    }
#endif
}

double local::DistortedPowerCorrelationFft::getPower(double k, double mu) const {
	if(mu < -1 || mu > 1) {
		throw RuntimeError("DistortedPowerCorrelationFft::getPower: expected -1 <= mu <= 1.");
	}
	return (*_power)(k)*(*_distortion)(k,mu);
}

double local::DistortedPowerCorrelationFft::getCorrelation(double r, double mu) const {
	if(mu < -1 || mu > 1) {
		throw RuntimeError("DistortedPowerCorrelationFft::getCorrelation: expected -1 <= mu <= 1.");
	}
	double rpar = r*mu;
	double rperp = r*std::sqrt(1-mu*mu);
	return (*_bicubicinterpolator)(rperp,rpar);
}

void local::DistortedPowerCorrelationFft::transform() {
#ifdef HAVE_LIBFFTW3
	// cleanup any previous plan
    if(_pimpl->data) FFTW(destroy_plan)(_pimpl->plan);
    // create a new plan for an in-place transform
    _pimpl->plan = FFTW(plan_dft_3d)(_nx,_ny,_nz,_pimpl->data,_pimpl->data,FFTW_BACKWARD,FFTW_ESTIMATE);
	// evaluate the power spectrum at each grid point (kx,ky,kz)
	for(int ix = 0; ix < _nx; ++ix){
		for(int iy = 0; iy < _ny; ++iy){
			for(int iz = 0; iz < _nz; ++iz){
				double ksq = _kxgrid[ix]*_kxgrid[ix] + _kygrid[iy]*_kygrid[iy] + _kzgrid[iz]*_kzgrid[iz];
            	double k = std::sqrt(ksq);
            	double mu = _kygrid[iy]/k;
            	std::size_t index(iz+_nz*(iy+_ny*ix));
            	_pimpl->data[index][0] = getPower(k,mu);
            	_pimpl->data[index][1] = 0;
            }
        }
    }
    // execute FFT to r space
	FFTW(execute)(_pimpl->plan);
	// extract the correlation function at grid points (rx,ry,0)
	for(int iy = 0; iy < _ny; ++iy) {
        for(int ix = 0; ix < _nx; ++ix) {
        	std::size_t ind(_nz*(iy+_ny*ix));
        	std::size_t ind2(ix+_nx*iy);
            _xi[ind2] = (_pimpl->data[ind][0]/_norm);
        }
    }
    // create the bicubic interpolator
	_bicubicinterpolator = new likely::BiCubicInterpolator(likely::BiCubicInterpolator::DataPlane(_xi),_spacing,_nx,_ny);
#endif
}
