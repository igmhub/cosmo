// Created 19-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_ADAPTIVE_MULTIPOLE_TRANSFORM
#define COSMO_ADAPTIVE_MULTIPOLE_TRANSFORM

#include "MultipoleTransform.h"

#include "likely/function.h"

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
		// Initializes for the specified function by automatically determining a suitable veps.
		// The termination criteria provided in the constructor will be tighted by a factor
		// 1/margin so that the nominal criteria are more likely to be met with other
		// similar functions.
		void initialize(likely::GenericFunctionPtr f, double margin = 2);
		// Calculates the transform of the specified function using the veps determined
		// from the most recent call to initialize(). Returns true if the termination
		// criteria are met. Results are stored in the vector provided, which will be
		// resized if necessary.
		bool transform(likely::GenericFunctionPtr f, std::vector<double> &results) const;
	private:
		MultipoleTransform::Type _type;
		int _ell;
		std::vector<double> _vpoints;
		double _relerr, _abserr, _abspow, _vmin, _vmax;
		bool _initialized;
	}; // AdaptiveMultipoleTransform
} // cosmo

#endif // COSMO_ADAPTIVE_MULTIPOLE_TRANSFORM
