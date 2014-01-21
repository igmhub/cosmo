// Created 20-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_DISTORTED_POWER_CORRELATION
#define COSMO_DISTORTED_POWER_CORRELATION

#include "cosmo/types.h"

namespace cosmo {
	class DistortedPowerCorrelation {
	// Represents the 3D correlation function corresponding to an isotropic power
	// spectrum that is distorted by a multiplicative function of (r,mu).
	public:
		DistortedPowerCorrelation(TabulatedPowerCPtr power, RMuFunctionCPtr distortion,
			int ellMax, double rmin, double rmax);
		virtual ~DistortedPowerCorrelation();
	private:
		TabulatedPowerCPtr _power;
		RMuFunctionCPtr _distortion;
		int _ellMax;
		double _rmin,_rmax;
	}; // DistortedPowerCorrelation
} // cosmo

#endif // COSMO_DISTORTED_POWER_CORRELATION
