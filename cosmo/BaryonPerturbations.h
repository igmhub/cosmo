// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_BARYON_PERTURBATIONS
#define COSMO_BARYON_PERTURBATIONS

#include "cosmo/types.h"

namespace cosmo {
    // Calculates baryon perturbations to a homogenous universe.
	class BaryonPerturbations {
	public:
		BaryonPerturbations(AbsHomogeneousUniversePtr homogeneous, double baryonFraction,
		    double hubbleConstant, double cmbTemperature);
		virtual ~BaryonPerturbations();
	private:
        AbsHomogeneousUniversePtr _homogeneous;
        double _baryonFraction, _hubbleConstant, _cmbTemperature;
	}; // BaryonPerturbations
} // cosmo

#endif // COSMO_BARYON_PERTURBATIONS
