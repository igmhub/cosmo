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
	// Allocate data array.
	_pimpl->data = (FFTW(complex)*) FFTW(malloc)(sizeof(FFTW(complex)) * nx*ny*nz);
#else
    throw RuntimeError("DistortedPowerCorrelationFft: package not built with FFTW3.");
#endif
	double twopi(8*std::atan(1));
    // Initialize the k space grid that will be used for evaluating the power spectrum.
    // Note that the grid is centered on (0,0,0) and grid points wrap around in a specific order.
    double dkx(twopi/(nx*spacing)), dky(twopi/(ny*spacing)), dkz(twopi/(nz*spacing));
	_kxgrid.reserve(nx);
	_kygrid.reserve(ny);
	_kzgrid.reserve(nz);
	for(int ix = 0; ix < nx; ++ix){
        double kx = (ix > nx/2 ? ix-nx : ix)*dkx;
        _kxgrid.push_back(kx);
    }
    for(int iy = 0; iy < ny; ++iy){
        double ky = (iy > ny/2 ? iy-ny : iy)*dky;
        _kygrid.push_back(ky);
    }
    for(int iz = 0; iz < nz; ++iz){
        double kz = (iz > nz/2 ? iz-nz : iz)*dkz;
        _kzgrid.push_back(kz);
    }
    // Calculate normalization.
    _norm = nx*ny*nz*spacing*spacing*spacing;
    // Initialize the array that will be used for bicubic interpolation.
    _xi.reset(new double[(nx/2+1)*(ny/2+1)]);
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
	// Clean up any previous plan.
    if(_pimpl->data) FFTW(destroy_plan)(_pimpl->plan);
    // Create a new plan for an in-place transform.
    _pimpl->plan = FFTW(plan_dft_3d)(_nx,_ny,_nz,_pimpl->data,_pimpl->data,FFTW_BACKWARD,FFTW_ESTIMATE);
	// Evaluate the power spectrum at each grid point (kx,ky,kz).
	for(int ix = 0; ix < _nx; ++ix){
		for(int iy = 0; iy < _ny; ++iy){
			for(int iz = 0; iz < _nz; ++iz){
				std::size_t index(iz+_nz*(iy+_ny*ix));
				double ksq = _kxgrid[ix]*_kxgrid[ix] + _kygrid[iy]*_kygrid[iy] + _kzgrid[iz]*_kzgrid[iz];
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
            	//if(iz==0 && ix==0) std::cout << k << " " << _pimpl->data[index][0] << std::endl;
            }
        }
    }
    // Execute FFT to r space.
	FFTW(execute)(_pimpl->plan);
	// Extract the correlation function at grid points (rx,ry,0).
	//double rvalue = -_spacing;
	for(int iy = 0; iy < _ny/2+1; ++iy) {
        for(int ix = 0; ix < _nx/2+1; ++ix) {
        	std::size_t ind(_nz*(iy+_ny*ix));
        	std::size_t ind2(ix+(_nx/2+1)*iy);
        	_xi[ind2] = _pimpl->data[ind][0]/_norm;
        	//if(iy==0) {
        	//	rvalue += _spacing;
        	//	std::cout << rvalue << " " << _pimpl->data[ind][0]/_norm << std::endl;
        	//}
        }
    }
    // Create the bicubic interpolator.
	_bicubicinterpolator = new likely::BiCubicInterpolator(likely::BiCubicInterpolator::DataPlane(_xi),_spacing,_nx/2+1,_ny/2+1);
#endif
}
