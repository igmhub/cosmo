// Created 09-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_POWER_SPECTRUM_CORRELATION_FUNCTION
#define COSMO_POWER_SPECTRUM_CORRELATION_FUNCTION

#include "cosmo/types.h"
#include "likely/types.h"

namespace cosmo {
    // Represents the l = 0,2,4 correlation functions of an isotropic 3D power
    // spectrum k^3/(2pi^2) P(k).
	class PowerSpectrumCorrelationFunction {
	public:
	    // Creates a correlation function from the power spectrum provided that is valid
	    // from rmin-rmax (in Mpc/h) and which uses nr logarithmically-spaced interpolation
	    // points over this range.
		PowerSpectrumCorrelationFunction(PowerSpectrumPtr powerSpectrum,
		    double rmin, double rmax, int nr = 1024);
		virtual ~PowerSpectrumCorrelationFunction();
		// Returns the correlation function evaluated at the specified radius in Mpc/h.
        double operator()(double rMpch) const;
	private:
        PowerSpectrumPtr _powerSpectrum; // evaluates k^3/(2pi^2) P(k)
        double _rmin, _rmax;
        int _nr;
        mutable likely::InterpolatorPtr _interpolator;
        mutable double _radius;
        // Integrand for 0 <= k < pi/r which is possibly singular and uses a series
        // expansion of sin(kr)/(kr) for small kr.
        double _integrand1(double kval) const;
        // Integrand for pi/r < k with the oscillatory sin(kr) part factored out.
        double _integrand2(double kval) const;
	}; // PowerSpectrumCorrelationFunction
} // cosmo

#endif // COSMO_POWER_SPECTRUM_CORRELATION_FUNCTION
