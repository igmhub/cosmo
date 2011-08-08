// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_HOMOGENEOUS_UNIVERSE_CALCULATOR
#define COSMO_HOMOGENEOUS_UNIVERSE_CALCULATOR

#include "cosmo/AbsHomogeneousUniverse.h"

#include "likely/types.h"

namespace cosmo {
    // Calculates the properties of a homogeneous and isotropic universe numerically
    // based on its Hubble function H(z)/H(0) and present-day curvature 1 - Omega(0).
	class HomogeneousUniverseCalculator : public AbsHomogeneousUniverse {
	public:
	    // Creates a new universe calculator that is valid for redshifts [0,zmax].
	    // Results are interpolated on a grid of nz redshifts spanning this range.
	    // Uses the specified target absolute accuracy for the interpolated points.
		HomogeneousUniverseCalculator(double zmax, int nz, double epsAbs);
		virtual ~HomogeneousUniverseCalculator();
		// Returns the present-day curvature defined as 1 - Omega(0).
        virtual double getCurvature() const = 0;
		// Returns the normalized Hubble function value H(z)/H(0) at the specified
		// redshift z >= 0.
        virtual double getHubbleFunction(double z) const = 0;
        // Returns the comoving line of sight distance in Mpc/h to an emitter with
        // the specified redshift z >= 0.
        virtual double getLineOfSightComovingDistance(double z) const;
        // Returns the comoving transverse distance scale in Mpc/h/rad between two
        // emitters at the same specified redshift z >= 0. Multiply this value by
        // the observed separation angle (rad) to obtain a physical distance in Mpc/h.
        virtual double getTransverseComovingScale(double z) const;
        // Returns the growth function D1(z) for small-scale perturbations in the absence
        // of neutrino free streaming. The result does not include the usual normalization
        // factor of 5/2*OmegaMatter (since we are not requiring that subclasses specify
        // a value of OmegaMatter and don't need one to do most of the calculation).
        virtual double getGrowthFunction(double z) const;
	private:
        double _zmax, _epsAbs;
        int _nz;
        mutable double _curvatureScale;
        mutable likely::InterpolatorPtr _lineOfSightInterpolator, _growthInterpolator;
        double _lineOfSightIntegrand(double z) const;
        double _growthIntegrand(double z) const;
	}; // HomogeneousUniverseCalculator
} // cosmo

#endif // COSMO_HOMOGENEOUS_UNIVERSE_CALCULATOR
