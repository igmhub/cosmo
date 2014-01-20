// Created 19-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/AdaptiveMultipoleTransform.h"
#include "cosmo/RuntimeError.h"

#include "likely/Interpolator.h"

#include "boost/foreach.hpp"

#include <cmath>
#include <iostream>

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

void local::AdaptiveMultipoleTransform::_evaluate(likely::GenericFunctionPtr f,
MultipoleTransformCPtr transform,std::vector<double> &result) {
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
		result[i] = interpolator(_vpoints[i]);
	}
}

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
		_mtGood.reset(new MultipoleTransform(_type, _ell, _vmin, _vmax, 2*vepsGuess,
			strategy, minSamplesPerCycle, minSamplesPerDecade, interpolationPadding));
		_evaluate(f,_mtGood,_resultsGood);
		_mtBetter.reset(new MultipoleTransform(_type, _ell, _vmin, _vmax, vepsGuess,
			strategy, minSamplesPerCycle, minSamplesPerDecade, interpolationPadding));
		_evaluate(f,_mtBetter,_resultsBetter);
		_veps = vepsGuess;
	}
	int tries(0), npoints(_vpoints.size());
	while(1) {
		std::cout << "try " << tries << ", veps = " << _veps << std::endl;
		// Give up after 100 tries
		if(++tries > 100) throw RuntimeError("AdaptiveMultipoleTransform: initialized failed.");
		// Check our termination criteria
		bool done(true);
		for(int i = 0; i < npoints; ++i) {
			double v(_vpoints[i]),fe(_resultsGood[i]),f2e(_resultsBetter[i]);
			double df = std::fabs(fe - f2e);
			std::cout << i << ' ' << v << ' ' << ' ' << fe << ' ' << f2e << ' ' << df << ' ' << df/std::fabs(f2e) << std::endl;
			if(df > _abserr*std::pow(v,_abspow)/margin && df > _relerr*std::fabs(f2e)/margin) {
				done = false;
				break;
			}
		}
		if(done) return _veps;
		// reduce veps by half and try again
		_mtGood = _mtBetter;
		_resultsGood.swap(_resultsBetter);
		_veps /= 2;
		_mtBetter.reset(new MultipoleTransform(_type, _ell, _vmin, _vmax, _veps,
			strategy, minSamplesPerCycle, minSamplesPerDecade, interpolationPadding));
		_evaluate(f,_mtBetter,_resultsBetter);
	}
}

bool local::AdaptiveMultipoleTransform::transform(
likely::GenericFunctionPtr f, std::vector<double> &result) const {
	if(!_mtGood || !_mtBetter) {
		throw RuntimeError("AdaptiveMultipoleTransform: must initialize before transforming.");
	}
	return true;
}
