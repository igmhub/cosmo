// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_RSD_POWER_SPECTRUM_CORRELATION_FUNCTION
#define COSMO_RSD_POWER_SPECTRUM_CORRELATION_FUNCTION

#include "cosmo/types.h"

namespace cosmo {
    // Represents the correlation functions of an isotropic 3D power
    // spectrum, input as k^3/(2pi^2) P(k), in redshift space (r,mu).
    // Anisotropic distortions are calculated in linear perturbation
    // theory, with the distant observer assumption. The original
    // reference for this is Kaiser, MNRAS 227, 1 (1987). See also
    // http://astro.berkeley.edu/~mwhite/teachdir/mini_red_mjw.pdf
	class RsdPowerSpectrumCorrelationFunction {
	public:
	    // Creates a correlation function from the power
	    // spectrum provided that is valid from rmin-rmax (in Mpc/h) and which uses
	    // nr logarithmically-spaced interpolation points over this range.
		RsdPowerSpectrumCorrelationFunction(PowerSpectrumPtr powerSpectrum,
		    double rmin, double rmax, int nr = 1024);
		virtual ~RsdPowerSpectrumCorrelationFunction();
	private:
	}; // RsdPowerSpectrumCorrelationFunction
} // cosmo

#endif // COSMO_RSD_POWER_SPECTRUM_CORRELATION_FUNCTION
