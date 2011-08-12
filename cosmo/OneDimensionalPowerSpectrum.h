// Created 10-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_ONE_DIMENSIONAL_POWER_SPECTRUM
#define COSMO_ONE_DIMENSIONAL_POWER_SPECTRUM

#include "cosmo/types.h"
#include "likely/types.h"

namespace cosmo {
    // Represents the one-dimensional projection of an isotropic 3D power spectrum.
	class OneDimensionalPowerSpectrum {
	public:
	    // Calculates the power spectrum of fluctuations along 1D lines corresponding
	    // to the specified 3D isotropic power spectrum. The resulting function is
	    // valid for the specified range of wavenumbers [kmin,kmax] and uses
	    // interpolation on nk logarithmically spaced points over this interval.
	    // The radius parameter specifies the thickness of the 1D to assume. A radius
	    // of zero corresponds to a mathematical zero-thickness line. A radius > 0 is
	    // interpreted as a hard cylindrical volume. A radius < 0 is interpreted as a
	    // soft Gaussian edge with sigma = -radius. Wavenumbers are in 1/(Mpc/h) and
	    // the radius is in Mpc/h.
		OneDimensionalPowerSpectrum(PowerSpectrumPtr powerSpectrum, double radius,
		    double kmin, double kmax, int nk = 1024);
		virtual ~OneDimensionalPowerSpectrum();
		// Returns the value of the one-dimensional power spectrum (kz/pi)P1(kz) for
		// the specified wavenumber in 1/(Mpc/h).
        double operator()(double kMpch) const;
	private:
        PowerSpectrumPtr _powerSpectrum;
        double _radius, _kmin, _kmax;
        int _nk;
        mutable likely::InterpolatorPtr _interpolator;
        mutable double _kz2;
        double _integrand(double kval) const;
	}; // OneDimensionalPowerSpectrum
} // cosmo

#endif // COSMO_ONE_DIMENSIONAL_POWER_SPECTRUM
