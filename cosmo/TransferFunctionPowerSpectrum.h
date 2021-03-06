// Created 09-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_TRANSFER_FUNCTION_POWER_SPECTRUM
#define COSMO_TRANSFER_FUNCTION_POWER_SPECTRUM

#include "cosmo/types.h"

#include "likely/function.h"

#include "boost/function.hpp"
#include "boost/smart_ptr.hpp"

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
	
	// The scoped global functions below provide generic power spectrum utilities...
	
	// Returns the RMS amplitude of fluctuations inside a sphere of the specified
	// radius in Mpc/h for the specified power spectrum function. Use gaussian = true
	// for a Gaussian window function with the specified radius. Otherwise, a step
	// function (top-hat) window function is used.
    double getRmsAmplitude(PowerSpectrumPtr powerSpectrum, double rMpch,
        bool gaussian = false);
	
	// Evaluates the Legendre polynomial for even ell up to 12 and returns 0 for any other ell.
    double legendreP(int ell, double mu);

    // Returns the specified multipole projection of the function provided, calculated
    // using numerical integration over 0 < mu < 1. Only even 0 <= ell <= 12 are implemented.
    double getMultipole(likely::GenericFunctionPtr fOfMuPtr, int ell, double epsAbs = 1e-6, double epsRel = 1e-6);
	
} // cosmo

#endif // COSMO_TRANSFER_FUNCTION_POWER_SPECTRUM
