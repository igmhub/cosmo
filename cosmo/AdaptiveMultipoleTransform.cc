// Created 19-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/AdaptiveMultipoleTransform.h"
#include "cosmo/RuntimeError.h"

#include "likely/Interpolator.h"

#include "boost/foreach.hpp"

#include <cmath>
#include <algorithm>

namespace local = cosmo;

local::AdaptiveMultipoleTransform::AdaptiveMultipoleTransform(MultipoleTransform::Type type,
int ell, double scale, std::vector<double>const &vpoints,
double relerr, double abserr, double abspow)
: _type(type), _ell(ell), _scale(scale), _vpoints(vpoints),
_relerr(relerr), _abserr(abserr), _abspow(abspow), _veps(0)
{
	// Input parameter validation
	if(_type != MultipoleTransform::SphericalBessel && _type != MultipoleTransform::Hankel) {
		throw RuntimeError("AdaptiveMultipoleTransform: invalid type.");
	}
	if(_ell < 0) {
		throw RuntimeError("AdaptiveMultipoleTransform: expected ell >= 0.");
	}
	if(_scale == 0) {
		throw RuntimeError("AdaptiveMultipoleTransform: expected scale != 0.");
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

void local::AdaptiveMultipoleTransform::_evaluate(likely::GenericFunctionPtr f,
MultipoleTransformCPtr transform,std::vector<double> &result) const {
	// Look up this transforms grids
	std::vector<double> const &ugrid = transform->getUGrid(), &vgrid = transform->getVGrid();
	// Prepare a grid of tabulated f(u) values
	std::vector<double> fgrid;
	fgrid.reserve(ugrid.size());
	BOOST_FOREACH(double u, ugrid) {
		fgrid.push_back((*f)(u));
	}
	// Calculate the corresponding grid of transform[f](v) values
	std::vector<double> ftgrid;
	(*transform).transform(fgrid,ftgrid);
	// Interpolate transform[f](v) to _vpoints
	likely::Interpolator interpolator(vgrid,ftgrid,"cspline");
	int npoints(_vpoints.size());
	if(result.size() != npoints) std::vector<double>(npoints).swap(result);
	for(int i = 0; i < npoints; ++i) {
		result[i] = _scale*interpolator(_vpoints[i]);
	}
}

bool local::AdaptiveMultipoleTransform::_isTerminated(double margin) const {
	for(int i = 0; i < _vpoints.size(); ++i) {
		double v(_vpoints[i]),f2e(_resultsGood[i]),fe(_resultsBetter[i]);
		double df = std::fabs(fe - f2e);
		if(df > _abserr*std::pow(v,_abspow)/margin && df > _relerr*std::fabs(fe)/margin) {
			return false;
		}
	}
	return true;
}

void local::AdaptiveMultipoleTransform::_saveResult(std::vector<double> &result) const {
	int npoints(_vpoints.size());
	if(result.size() != npoints) {
		// Replace results with a vector of the required size
		std::vector<double>(npoints).swap(result);
	}
	// Swap the contents of result and our internal result vector. This just swaps
	// pointers so is fast, and leaves our internal result vector with undefined contents.
	_resultsBetter.swap(result);
}

double local::AdaptiveMultipoleTransform::initialize(
likely::GenericFunctionPtr f, std::vector<double> &result,
int minSamplesPerDecade, double margin, double vepsMax, double vepsMin, bool optimize) {
	if(margin < 1) {
		throw RuntimeError("AdaptiveMultipoleTransform: expected margin >= 1.");
	}
	MultipoleTransform::Strategy strategy(MultipoleTransform::EstimatePlan);
	int minSamplesPerCycle(2),interpolationPadding(3);
	// Create our first pair of transformers, if necessary
	if(!_mtGood || !_mtBetter) {
		if(vepsMax <= 0) {
			throw RuntimeError("AdaptiveMultipoleTransform: expected vepsGuess > 0.");
		}
		// Find an initial veps starting from vepsMax.
		_veps = vepsMax;
		while(_veps > vepsMin) {
			// Initialize without any min samples per decade, so we can see what it
			// would be for this trial veps
			int noMinSamplesPerDecade(0);
			_mtBetter.reset(new MultipoleTransform(_type, _ell, _vmin, _vmax, _veps,
				strategy, minSamplesPerCycle, noMinSamplesPerDecade, interpolationPadding));
			// Is this veps small enough to meet our samples/decade requirement?
			if(_mtBetter->getSamplesPerDecade() >= minSamplesPerDecade) break;
			// Otherwise, try a smaller veps
			_veps /= 2;
		}
		// Calculate the corresponding prediction
		_evaluate(f,_mtBetter,_resultsBetter);
		// Initialize a "good" transformer with veps that is 2x larger
		_mtGood.reset(new MultipoleTransform(_type, _ell, _vmin, _vmax, 2*_veps,
			strategy, minSamplesPerCycle, minSamplesPerDecade, interpolationPadding));
		_evaluate(f,_mtGood,_resultsGood);
	}
	while(_veps > vepsMin) {
		// Check our termination criteria
		if(_isTerminated(margin)) {
			_saveResult(result);
			if(optimize) {
				// Recreate transform objects using the MeasurePlan strategy
				strategy = MultipoleTransform::MeasurePlan;
				_mtGood.reset(new MultipoleTransform(_type, _ell, _vmin, _vmax, 2*_veps,
					strategy, minSamplesPerCycle, minSamplesPerDecade, interpolationPadding));
				_mtBetter.reset(new MultipoleTransform(_type, _ell, _vmin, _vmax, _veps,
					strategy, minSamplesPerCycle, minSamplesPerDecade, interpolationPadding));
			}
			return _veps;
		}
		// reduce veps by half and try again
		_veps /= 2;
		if(_veps < vepsMin) {
			throw RuntimeError("AdaptiveMultipoleTransform: reached vepsMin without convergence.");
		}
		_mtGood = _mtBetter;
		_resultsGood.swap(_resultsBetter);
		_mtBetter.reset(new MultipoleTransform(_type, _ell, _vmin, _vmax, _veps,
			strategy, minSamplesPerCycle, minSamplesPerDecade, interpolationPadding));
		_evaluate(f,_mtBetter,_resultsBetter);
	}
}

bool local::AdaptiveMultipoleTransform::transform(
likely::GenericFunctionPtr f, std::vector<double> &result, bool bypassTerminationTest) const {
	if(!_mtGood || !_mtBetter) {
		throw RuntimeError("AdaptiveMultipoleTransform: must initialize before transforming.");
	}
	_evaluate(f,_mtBetter,_resultsBetter);
	bool accurate(true);
	if(!bypassTerminationTest) {
		_evaluate(f,_mtGood,_resultsGood);
		accurate = _isTerminated();
	}
	_saveResult(result);
	return accurate;
}

double local::AdaptiveMultipoleTransform::getUMin() const {
	if(!_mtBetter) {
		throw RuntimeError("AdaptiveMultipoleTransform::getUMin: must initialize first.");
	}
	// ugrid is decreasing, so last element is the min value
	return _mtBetter->getUGrid().back();
}

double local::AdaptiveMultipoleTransform::getUMax() const {
	if(!_mtBetter) {
		throw RuntimeError("AdaptiveMultipoleTransform::getUMax: must initialize first.");
	}
	// ugrid is decreasing, so first element is the min value
	return _mtBetter->getUGrid().front();
}

int local::AdaptiveMultipoleTransform::getNU() const {
	if(!_mtBetter) {
		throw RuntimeError("AdaptiveMultipoleTransform::getNU: must initialize first.");
	}
	return _mtBetter->getUGrid().size();
}

double local::AdaptiveMultipoleTransform::getUSamplesPerDecade() const {
	return getNU()/std::log10(getUMax()/getUMin());
}