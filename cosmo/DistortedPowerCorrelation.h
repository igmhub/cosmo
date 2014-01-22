// Created 20-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_DISTORTED_POWER_CORRELATION
#define COSMO_DISTORTED_POWER_CORRELATION

#include "cosmo/types.h"
#include "likely/function.h"

#include "boost/smart_ptr.hpp"

#include <vector>

namespace cosmo {
	class AdaptiveMultipoleTransform;
	class DistortedPowerCorrelation {
	// Represents the 3D correlation function corresponding to an isotropic power
	// spectrum that is distorted by a multiplicative function of (r,mu).
	public:
		// Creates a new distorted power correlation function using the specified
		// isotropic power P(k) and distortion function D(k,mu). The resulting
		// correlation function will be valid over [rmin,rmax] and estimated
		// using multipoles up to ellMax that are interpolated using nr equally-spaced
		// points in r. If symmetric is true, then only even multipoles are used.
		// The desired accuracy is specified by relerr, abserr, and abspow, such that
		// the difference between the true and estimated xi(r,mu) satisfies:
		// |true-est| < max(abserr*r^abspow,true*true)
		DistortedPowerCorrelation(likely::GenericFunctionPtr power, RMuFunctionCPtr distortion,
			double rmin, double rmax, int nr, int ellMax, bool symmetric = true,
			double relerr = 1e-2, double abserr = 1e-3, double abspow = 0);
		virtual ~DistortedPowerCorrelation();
		// Returns the value of P(k,mu) = P(k)*D(k,mu)
		double getPower(double k, double mu) const;
		// Returns the specified multipole of P(k,mu) evaluated at k
		double getPowerMultipole(double k, int ell) const;
		// Initializes our multipole estimates and correlation transforms.
		void initialize();
		// Transforms the k-space power multipoles to r space. Returns true if the termination
		// criteria are met, unless bypassTerminationTest is true (in which case we
		// always return true and transforms will be somewhat faster).
		bool transform(bool bypassTerminationTest = false) const;
	private:
		likely::GenericFunctionPtr _power;
		RMuFunctionCPtr _distortion;
		double _relerr,_abserr,_abspow;
		int _ellMax;
		bool _symmetric;
		std::vector<double> _rgrid;
		mutable std::vector<std::vector<double> > _xiMoments;
		typedef boost::shared_ptr<AdaptiveMultipoleTransform> AdaptiveMultipoleTransformPtr;
		std::vector<AdaptiveMultipoleTransformPtr> _transformer;
	}; // DistortedPowerCorrelation
} // cosmo

#endif // COSMO_DISTORTED_POWER_CORRELATION
