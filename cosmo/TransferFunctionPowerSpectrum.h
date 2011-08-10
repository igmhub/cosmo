// Created 09-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_TRANSFER_FUNCTION_POWER_SPECTRUM
#define COSMO_TRANSFER_FUNCTION_POWER_SPECTRUM

#include "cosmo/types.h"

namespace cosmo {
    // Represents an isotropic power spectrum of 3D inhomogeneities based on a model of
    // primordial fluctuations and a transfer function.
	class TransferFunctionPowerSpectrum {
	public:
		TransferFunctionPowerSpectrum(TransferFunctionPtr transferFunction,
		    double spectralIndex = 1, double deltaH = 1);
		virtual ~TransferFunctionPowerSpectrum();
		// Returns the value of k^3/(2pi^2) P(k).
        double operator()(double kMpch) const;
	private:
        TransferFunctionPtr _transferFunction;
        double _spectralIndex, _deltaHSq;
	}; // TransferFunctionPowerSpectrum
	
	// Returns the RMS amplitude of fluctuations inside a sphere of the specified
	// radius in Mpc/h for the specified power spectrum function. Use gaussian = true
	// for a Gaussian window function with the specified radius. Otherwise, a step
	// function (top-hat) window function is used.
    double getRmsAmplitude(PowerSpectrumPtr powerSpectrum, double rMpch,
        bool gaussian = false);
	
} // cosmo

#endif // COSMO_TRANSFER_FUNCTION_POWER_SPECTRUM
