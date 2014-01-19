// Created 19-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_ADAPTIVE_MULTIPOLE_TRANSFORM
#define COSMO_ADAPTIVE_MULTIPOLE_TRANSFORM

#include "MultipoleTransform.h"

#include <vector>

namespace cosmo {
	class AdaptiveMultipoleTransform {
	// Uses the MultipoleTransform class to calculate transforms but replaces the
	// veps numerical control parameter with accuracy criteria that are used to
	// automatically set veps. Also takes a function pointer as input, instead of
	// requiring the user to tabulate function values.
	public:
		// Creates a new transformer of the specified type and multipole. Subsequent
		// transforms will be provided at the specified vpoints, which will also be
		// used to adaptively monitor numerical errors. The numerical termination
		// criteria is that |f(2*veps) - f(veps)| < max(abserr*v^abspow,relerr*|f(veps)|)
		// for each point v in vpoints.
		AdaptiveMultipoleTransform(MultipoleTransform::Type type, int ell,
			std::vector<double> const &vpoints, double relerr, double abserr, double abspow = 0);
		virtual ~AdaptiveMultipoleTransform();
	private:
		MultipoleTransform::Type _type;
		int _ell;
		std::vector<double> _vpoints;
		double _relerr, _abserr, _abspow, _vmin, _vmax;
	}; // AdaptiveMultipoleTransform
} // cosmo

#endif // COSMO_ADAPTIVE_MULTIPOLE_TRANSFORM
