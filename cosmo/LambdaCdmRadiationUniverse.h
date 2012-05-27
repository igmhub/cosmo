// Created 26-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_LAMBDA_CDM_RADIATION_UNIVERSE
#define COSMO_LAMBDA_CDM_RADIATION_UNIVERSE

#include "cosmo/HomogeneousUniverseCalculator.h"

namespace cosmo {
	class LambdaCdmRadiationUniverse : public HomogeneousUniverseCalculator {
	public:
		LambdaCdmRadiationUniverse(double OmegaMatter, double OmegaK = 0,
		    double hValue = 0.7, double Tcmb = 2.725, double Nnu = 3.046,
		    double zmax = 10, int nz = 1000, double epsAbs = 1e-8);
		virtual ~LambdaCdmRadiationUniverse();
		// Returns the present-day curvature defined as 1 - Omega(0).
        virtual double getCurvature() const;
		// Returns the normalized Hubble function value H(z)/H(0) at the specified
		// redshift z >= 0.
        virtual double getHubbleFunction(double z) const;
	private:
        double _OmegaMatter, _OmegaK, _OmegaRadiation, _OmegaLambda;
	}; // LambdaCdmRadiationUniverse
} // cosmo

#endif // COSMO_LAMBDA_CDM_RADIATION_UNIVERSE
