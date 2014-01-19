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

double local::AdaptiveMultipoleTransform::initialize(
likely::GenericFunctionPtr f, int minSamplesPerDecade, double margin, double vepsGuess) {
	if(margin < 1) {
		throw RuntimeError("AdaptiveMultipoleTransform: expected margin >= 1.");
	}
	MultipoleTransform::Strategy strategy(MultipoleTransform::EstimatePlan);
	int minSamplesPerCycle(2),interpolationPadding(3);
	// Create our first pair of transformers, if necessary
	if(!_mtGood || !_mtBetter) {
		if(vepsGuess <= 0) {
			throw RuntimeError("AdaptiveMultipoleTransform: expected vepsGuess > 0.");
		}
		_mtBetter.reset(new MultipoleTransform(_type, _ell, _vmin, _vmax, vepsGuess,
			strategy, minSamplesPerCycle, minSamplesPerDecade, interpolationPadding));
		_mtGood.reset(new MultipoleTransform(_type, _ell, _vmin, _vmax, 2*vepsGuess,
			strategy, minSamplesPerCycle, minSamplesPerDecade, interpolationPadding));
		_veps = vepsGuess;
	}
	return _veps;
}

bool local::AdaptiveMultipoleTransform::transform(
likely::GenericFunctionPtr f, std::vector<double> &results) const {
	if(!_mtGood || !_mtBetter) {
		throw RuntimeError("AdaptiveMultipoleTransform: must initialize before transforming.");
	}
	return true;
}
