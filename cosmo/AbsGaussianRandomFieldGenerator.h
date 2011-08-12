// Created 12-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_ABS_GAUSSIAN_RANDOM_FIELD_GENERATOR
#define COSMO_ABS_GAUSSIAN_RANDOM_FIELD_GENERATOR

#include "cosmo/types.h"

namespace cosmo {
    // Represents an abstract generator of 3D Gaussian random fields as realizations of
    // a 3D isotropic power spectrum on a uniform rectangular grid.
	class AbsGaussianRandomFieldGenerator {
	public:
	    // Creates a new generator with the specified grid spacing in Mpc/h and the
	    // specified dimensions (nx,ny,nz). No realization is generated until the
	    // generate method is called.
		AbsGaussianRandomFieldGenerator(PowerSpectrumPtr powerSpectrum, double spacing,
		    int nx, int ny, int nz);
		virtual ~AbsGaussianRandomFieldGenerator();
		// Accessors for constructor parameters.
        double getSpacing() const;
        int getNx() const;
        int getNy() const;
        int getNz() const;
        // Generates a new realization using the specified random seed and stores the
        // results internally. Use the getField() method to access generated values.
        void generate(int seed);
        // Returns the most recent generated value at the specified grid point. Throws
        // a RuntimeError if generate() has never been called.
        double getField(int x, int y, int z) const;
        // Returns the number of times that generate() has been called.
        int getGenerateCount() const;
	private:
        PowerSpectrumPtr _powerSpectrum;
        double _spacing;
        int _nx,_ny,_nz, _generateCount;
        // The generate() method calls this virtual method.
        virtual void _generate(int seed) = 0;
        // The getField method calls this virtual method after checking that data is
        // available and that (x,y,z) are valid values.
        virtual double _getFieldUnchecked(int x, int y, int z) const = 0;
	}; // AbsGaussianRandomFieldGenerator
	
    inline double AbsGaussianRandomFieldGenerator::getSpacing() const { return _spacing; }

    inline int AbsGaussianRandomFieldGenerator::getNx() const { return _nx; }
    inline int AbsGaussianRandomFieldGenerator::getNy() const { return _ny; }
    inline int AbsGaussianRandomFieldGenerator::getNz() const { return _nz; }

    inline int AbsGaussianRandomFieldGenerator::getGenerateCount() const {
        return _generateCount;
    }
    inline void AbsGaussianRandomFieldGenerator::generate(int seed) {
        _generateCount++;
        _generate(seed);
    }
	
} // cosmo

#endif // COSMO_ABS_GAUSSIAN_RANDOM_FIELD_GENERATOR
