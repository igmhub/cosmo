// Created 12-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_FFT_GAUSSIAN_RANDOM_FIELD_GENERATOR
#define COSMO_FFT_GAUSSIAN_RANDOM_FIELD_GENERATOR

#include "cosmo/AbsGaussianRandomFieldGenerator.h"

#include "boost/smart_ptr.hpp"

namespace cosmo {
    // Implements the abstract Gaussian random field generator interface using FFT.
	class FftGaussianRandomFieldGenerator : public AbsGaussianRandomFieldGenerator {
	public:
		FftGaussianRandomFieldGenerator(PowerSpectrumPtr powerSpectrum, double spacing,
		    int nx, int ny, int nz, likely::RandomPtr random = likely::RandomPtr());
		virtual ~FftGaussianRandomFieldGenerator();
        // Generates a new realization and stores the results internally. Use the getField()
        // method to access generated values.
        virtual void generate();
        // Returns the memory size in bytes required for this generator or zero if this
        // information is not available.
        virtual std::size_t getMemorySize() const;
	private:
        class Implementation;
        int _halfz;
        std::size_t _nbuf;
        boost::scoped_ptr<Implementation> _pimpl;
        // The getField method calls this after checking for invalid (x,y,z).
        virtual double _getFieldUnchecked(int x, int y, int z) const;
	}; // FftGaussianRandomFieldGenerator
} // cosmo

#endif // COSMO_FFT_GAUSSIAN_RANDOM_FIELD_GENERATOR
