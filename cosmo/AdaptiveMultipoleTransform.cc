// Created 19-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/AdaptiveMultipoleTransform.h"
#include "cosmo/RuntimeError.h"

namespace local = cosmo;

local::AdaptiveMultipoleTransform::AdaptiveMultipoleTransform(MultipoleTransform::Type type,
int ell, std::vector<double>const &vpoints, double relerr, double abserr, double abspow)
: _type(type), _ell(ell), _vpoints(vpoints), _relerr(relerr), _abserr(abserr), _abspow(abspow)
{
	// Input parameter validation
	if(_type != MultipoleTransform::SphericalBessel && _type != MultipoleTransform::Hankel) {
		throw RuntimeError("AdaptiveMultipoleTransform: invalid type.");
	}
	if(_ell < 0) {
		throw RuntimeError("AdaptiveMultipoleTransform: expected ell >= 0.");
	}
	if(_relerr <= 0 && _abserr <= 0) {
		throw RuntimeError("AdaptiveMultipoleTransform: invalid termination criteria.");
	}
	int nv(_vpoints.size());
	if(nv < 2) {
		throw RuntimeError("AdaptiveMultipoleTransform: expected at least 2 vpoints.");
	}
	for(int iv = 1; iv < nv; ++iv) {
		if(_vpoints[iv] <= _vpoints[iv-1]) {
			throw RuntimeError("AdaptiveMultipoleTransform: vpoints not increasing.");
		}
	}
	_vmin = _vpoints.front();
	_vmax = _vpoints.back();
	if(_vmin <= 0) {
		throw RuntimeError("AdaptiveMultipoleTransform: expected vmin > 0.");
	}
}

local::AdaptiveMultipoleTransform::~AdaptiveMultipoleTransform() { }
