// Created 20-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_DISTORTED_POWER_CORRELATION
#define COSMO_DISTORTED_POWER_CORRELATION

#include "cosmo/types.h"
#include "likely/types.h"
#include "likely/function.h"

#include "boost/smart_ptr.hpp"

#include <vector>
#include <iosfwd>

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
		// Returns the value of P(k,mu) = P(k)*D(k,mu). This is fast to evaluate and
		// does not require that initialize() be called first.
		double getPower(double k, double mu) const;
		// Returns the specified multipole of P(k,mu) evaluated at k. This method calculates
		// the relevant numerical 1D integral each time it is called, so is relatively slow,
		// but does not require that initialize() be called first.
		double getPowerMultipole(double k, int ell) const;
		// Returns the specified multipole of xi(r,mu) evaluated at r. This is relatively fast
		// to evaluate, only involving some interpolation, but requires that initialize()
		// be called first.
		double getCorrelationMultipole(double r, int ell) const;
		// Returns the correlation function xi(r,mu). This is relatively fast
		// to evaluate, only involving some interpolation, but requires that initialize()
		// be called first.
		double getCorrelation(double r, double mu) const;
		// Initializes our multipole estimates and correlation transforms and automatically
		// sets the relerr and abserr goals for each multipole based on their relative
		// contributions in [rmin,rmax], determined by sampling a nr-by-nmu grid. For other
		// options, see AdaptiveMultipoleTransform::initialize().
		void initialize(int nmu = 20, int minSamplesPerDecade= 40, double margin = 2,
			double vepsMax = 0.01, double vepsMin = 1e-6, bool optimize = false);
		// Tests if we have ever been initialized.
		bool isInitialized() const;
		// Transforms the k-space power multipoles to r space. Returns true if the termination
		// criteria are met, unless bypassTerminationTest is true (in which case we
		// always return true and transforms will be somewhat faster).
		bool transform(bool bypassTerminationTest = false) const;
		// Returns a shared const pointer to the specified transform.
		AdaptiveMultipoleTransformCPtr getTransform(int ell) const;
		// Fills the variables provided with the (r,mu) coordinates where the specified
		// multipole has the biggest relative contrbution:
		//
		//   rel = |f_ell(r)*L_ell(mu)|/|xi(r,mu)|
		//
		// to the estimated correlation function, as well as the value of rel at (r,mu).
		void getBiggestContribution(int ell, double &rbig, double &mubig, double &relbig) const;
		// Prints info about this object to the specified output stream.
		void printToStream(std::ostream &out) const;
	private:
		likely::GenericFunctionPtr _power;
		RMuFunctionCPtr _distortion;
		double _relerr,_abserr,_abspow;
		int _ellMax;
		bool _symmetric, _initialized;
		std::vector<double> _rgrid, _rbig, _mubig, _relbig;
		mutable std::vector<std::vector<double> > _xiMoments;
		mutable std::vector<likely::InterpolatorPtr> _interpolator;
		std::vector<AdaptiveMultipoleTransformPtr> _transformer;
	}; // DistortedPowerCorrelation

	inline bool DistortedPowerCorrelation::isInitialized() const { return _initialized; }

} // cosmo

#endif // COSMO_DISTORTED_POWER_CORRELATION
