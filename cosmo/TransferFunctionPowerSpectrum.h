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
		// Gets/sets the value of spectralIndex.
        double getSpectralIndex() const;
        void setSpectralIndex(double value);
        // Gets/sets the value of deltaH.
        double getDeltaH() const;
        void setDeltaH(double value);
		// Returns the value of k^3/(2pi^2) P(k).
        double operator()(double kMpch) const;
        // Sets our normalization (deltaH) to match the specified value of sigma, the
        // RMS amplitude of fluctuations within a radius of rMpch, and returns the
        // previous value. Use gaussian = true for a Gaussian window function with the
        // specified radius. Otherwise, a step function (top-hat) window function is
        // used. The default values of rMpch = 8 and gaussian = false correspond to the
        // usual definition of sigma8.
        double setSigma(double sigma, double rMpch = 8, bool gaussian = false);
	private:
        TransferFunctionPtr _transferFunction;
        double _spectralIndex, _deltaH, _deltaHSq;
	}; // TransferFunctionPowerSpectrum
	
    inline double TransferFunctionPowerSpectrum::getSpectralIndex() const { return _spectralIndex; }
    inline double TransferFunctionPowerSpectrum::getDeltaH() const { return _deltaH; }
	
	// Returns the RMS amplitude of fluctuations inside a sphere of the specified
	// radius in Mpc/h for the specified power spectrum function. Use gaussian = true
	// for a Gaussian window function with the specified radius. Otherwise, a step
	// function (top-hat) window function is used.
    double getRmsAmplitude(PowerSpectrumPtr powerSpectrum, double rMpch,
        bool gaussian = false);
	
} // cosmo

#endif // COSMO_TRANSFER_FUNCTION_POWER_SPECTRUM
