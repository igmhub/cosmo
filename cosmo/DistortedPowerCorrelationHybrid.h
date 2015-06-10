// Created 11-May-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#ifndef COSMO_DISTORTED_POWER_CORRELATION_HYBRID
#define COSMO_DISTORTED_POWER_CORRELATION_HYBRID

#include "cosmo/types.h"
#include "likely/types.h"
#include "likely/function.h"

#include "boost/smart_ptr.hpp"

#include <vector>
#include <iosfwd>

namespace likely{class BiCubicInterpolator;}
namespace cosmo {
	class DistortedPowerCorrelationHybrid {
	// Represents the 3D correlation function corresponding to an isotropic power
	// spectrum P(k) that is distorted by a multiplicative function D(k,mu_k).
	// This class is optimized for the case where D(k,mu_k) is allowed to
	// change (e.g., as its internal parameters are changed) but that after
	// each change, the correlation function xi(r,mu) needs to be evaluated
	// many times. The normal usage is:
	//
	//    - transform() each time D(k,mu_k) changes internally
	//      - call getCorrelation(r,mu) many times
	//
	// The transformation employs a two-step algorithm that first evaluates a series
	// of 1D Fourier transformations, followed by a series of 1D integrals.
	// Cartesian x axis represents the transverse direction, and cartesian y axis
	// represents the radial direction.
	public:
		// Creates a new distorted power correlation function using the specified
		// isotropic power P(k) and distortion function D(k,mu).
		DistortedPowerCorrelationHybrid(likely::GenericFunctionPtr power, KMuPkFunctionCPtr distortion,
		    double kxmin, double kxmax, int nx, double spacing, int ny, double rmax, double epsAbs = 1e-6);
		virtual ~DistortedPowerCorrelationHybrid();
		// Returns the value of P(k,mu) = P(k)*D(k,mu).
		double getPower(double k, double mu) const;
		// Returns the correlation function xi(r,mu).
		double getCorrelation(double r, double mu) const;
		// Returns the k-space transform ktf(ry,kx).
		double getKTransform(double ry, double kx) const;
		// Transforms the k-space power spectrum to r-space.
		void transform();
		// Performs a series of 1D Fourier transforms of k-space power spectrum.
		void ktransform();
		// Returns the memory size in bytes required for this transform or zero if this
        // information is not available.
        virtual std::size_t getMemorySize() const;
	private:
		class Implementation;
		boost::scoped_ptr<Implementation> _pimpl;
		likely::GenericFunctionPtr _power;
		KMuPkFunctionCPtr _distortion;
		std::vector<double> _kxgrid, _kygrid, _rgrid;
		boost::shared_array<double> _ktf, _xi;
		double _kxmin, _kxmax, _spacing, _rmax, _epsAbs, _twopi, _norm, _rx, _dkx;
		int _nx, _ny, _nr, _count;
		double _transverseIntegrand(double kx) const;
		likely::BiCubicInterpolator *_xiInterpolator, *_ktransformInterpolator;
		mutable likely::InterpolatorPtr _ktfInterpolator;
	}; // DistortedPowerCorrelationHybrid

} // cosmo

#endif // COSMO_DISTORTED_POWER_CORRELATION_HYBRID
