// Created 12-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_ABS_GAUSSIAN_RANDOM_FIELD_GENERATOR
#define COSMO_ABS_GAUSSIAN_RANDOM_FIELD_GENERATOR

#include "cosmo/types.h"
#include "likely/types.h"

#include <cstddef>

namespace cosmo {
    // Represents an abstract generator of 3D Gaussian random fields as realizations of
    // a 3D isotropic power spectrum on a uniform rectangular grid.
	class AbsGaussianRandomFieldGenerator {
	public:
	    // Creates a new generator with the specified grid spacing in Mpc/h and the
	    // specified dimensions (nx,ny,nz). No realization is generated until the
	    // first call to generate(). Uses the random number source provided, or else
	    // the default likely::Random::instance().
		AbsGaussianRandomFieldGenerator(PowerSpectrumPtr powerSpectrum, double spacing,
		    int nx, int ny, int nz, likely::RandomPtr random = likely::RandomPtr());
		virtual ~AbsGaussianRandomFieldGenerator();
		// Accessors for constructor parameters.
        double getSpacing() const;
        int getNx() const;
        int getNy() const;
        int getNz() const;
        // Generates a new realization and stores the results internally. Use the getField()
        // method to access generated values.
        virtual void generate() = 0;
        // Returns the most recent generated value at the specified grid point. Throws
        // a RuntimeError for invalid (x,y,z).
        double getField(int x, int y, int z) const;
        // Returns the memory size in bytes required for this generator or zero if this
        // information is not available.
        virtual std::size_t getMemorySize() const;
    protected:
        likely::RandomPtr getRandom();
	private:
        PowerSpectrumPtr _powerSpectrum;
        double _spacing;
        int _nx,_ny,_nz;
        likely::RandomPtr _random;
        // The getField method calls this virtual method after checking that data is
        // available and that (x,y,z) are valid values.
        virtual double _getFieldUnchecked(int x, int y, int z) const = 0;
	}; // AbsGaussianRandomFieldGenerator
	
    inline double AbsGaussianRandomFieldGenerator::getSpacing() const { return _spacing; }

    inline int AbsGaussianRandomFieldGenerator::getNx() const { return _nx; }
    inline int AbsGaussianRandomFieldGenerator::getNy() const { return _ny; }
    inline int AbsGaussianRandomFieldGenerator::getNz() const { return _nz; }

    inline likely::RandomPtr AbsGaussianRandomFieldGenerator::getRandom() { return _random; }
	
} // cosmo

#endif // COSMO_ABS_GAUSSIAN_RANDOM_FIELD_GENERATOR
