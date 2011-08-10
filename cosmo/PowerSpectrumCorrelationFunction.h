// Created 09-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_POWER_SPECTRUM_CORRELATION_FUNCTION
#define COSMO_POWER_SPECTRUM_CORRELATION_FUNCTION

#include "cosmo/types.h"

namespace cosmo {
    // Represents the isotropic correlation function of an isotropic 3D power spectrum.
	class PowerSpectrumCorrelationFunction {
	public:
	    // Creates a correlation function from the power spectrum provided that is valid
	    // from rmin-rmax (in Mpc/h) and which uses nr equally-spaced interpolation points
	    // over this range.
		PowerSpectrumCorrelationFunction(PowerSpectrumPtr powerSpectrum,
		    double rmin, double rmax, int nr = 8192);
		virtual ~PowerSpectrumCorrelationFunction();
		// Returns the correlation function evaluated at the specified radius in Mpc/h.
        double operator()(double rMpch) const;
	private:
        PowerSpectrumPtr _powerSpectrum;
        double _rmin, _rmax;
        int _nr;
	}; // PowerSpectrumCorrelationFunction
} // cosmo

#endif // COSMO_POWER_SPECTRUM_CORRELATION_FUNCTION
