// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_LAMBDA_CDM_UNIVERSE
#define COSMO_LAMBDA_CDM_UNIVERSE

#include "cosmo/HomogeneousUniverseCalculator.h"

namespace cosmo {
	class LambdaCdmUniverse : public HomogeneousUniverseCalculator {
	public:
		LambdaCdmUniverse(double OmegaLambda, double OmegaMatter,
		    double zmax = 10, int nz = 1000, double epsAbs = 1e-8);
		virtual ~LambdaCdmUniverse();
		// Returns the present-day curvature defined as 1 - Omega(0).
        virtual double getCurvature() const;
		// Returns the normalized Hubble function value H(z)/H(0) at the specified
		// redshift z >= 0.
        virtual double getHubbleFunction(double z) const;
	private:
        double _OmegaLambda, _OmegaMatter, _curvature;
	}; // LambdaCdmUniverse
} // cosmo

#endif // COSMO_LAMBDA_CDM_UNIVERSE
