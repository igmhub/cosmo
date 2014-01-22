// Created 19-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_ADAPTIVE_MULTIPOLE_TRANSFORM
#define COSMO_ADAPTIVE_MULTIPOLE_TRANSFORM

#include "cosmo/MultipoleTransform.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"

#include "likely/function.h"

#include "boost/smart_ptr.hpp"

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
		// similar functions. If this is the first time we have been initialized, uses
		// vepsMax as a starting point. Otherwise, the veps value from the initialization is
		// used as the starting point. The veps value is then successively halved until the
		// termination criteria are met or we hit the vepsMin limit (which throws a
		// RuntimeError). Results are stored in the vector provided, which will be
		// resized if necessary. If optimize is true, then we perform an additional step
		// of optimizing the FFTs that will be necessary for subsequent transforms. This
		// optimization step takes at least a few seconds so is only worth doing if
		// many transforms will be performed per initialization. Note that optimized
		// transforms will generally give different numerical results at the level of
		// roundoff errors. Returns the selected veps value.
		double initialize(likely::GenericFunctionPtr f, std::vector<double> &result,
			int minSamplesPerDecade= 40, double margin = 2,
			double vepsMax = 0.01, double vepsMin = 1e-6, bool optimize = false);
		// Calculates the transform of the specified function using the veps determined
		// from the most recent call to initialize(). Results are stored in the vector
		// provided, which will be resized if necessary. Returns true if the termination
		// criteria are met, unless bypassTerminationTest is true (in which case we
		// always return true and transforms will be somewhat faster).
		bool transform(likely::GenericFunctionPtr f, std::vector<double> &result,
			bool bypassTerminationTest = false) const;
	private:
		MultipoleTransform::Type _type;
		int _ell;
		std::vector<double> _vpoints;
		mutable std::vector<double> _resultsGood, _resultsBetter;
		double _relerr, _abserr, _abspow, _vmin, _vmax, _veps;
		typedef boost::shared_ptr<const MultipoleTransform> MultipoleTransformCPtr;
		MultipoleTransformCPtr _mtGood, _mtBetter;
		void _evaluate(likely::GenericFunctionPtr f,
			MultipoleTransformCPtr transform, std::vector<double> &result) const;
		bool _isTerminated(double margin = 1) const;
		void _saveResult(std::vector<double> &result) const;
	}; // AdaptiveMultipoleTransform
} // cosmo

#endif // COSMO_ADAPTIVE_MULTIPOLE_TRANSFORM
