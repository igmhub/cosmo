// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/BaryonPerturbations.h"
#include "cosmo/RuntimeError.h"
#include "cosmo/AbsHomogeneousUniverse.h"

namespace local = cosmo;

local::BaryonPerturbations::BaryonPerturbations(AbsHomogeneousUniversePtr homogeneous,
double baryonFraction, double hubbleConstant, double cmbTemperature)
: _homogeneous(homogeneous), _baryonFraction(baryonFraction), _hubbleConstant(hubbleConstant),
_cmbTemperature(cmbTemperature)
{
    if(!homogeneous) {
        throw RuntimeError("BaryonPerturbation: missing required homogeneous universe.");
    }
    if(baryonFraction < 0 || baryonFraction > 1) {
        throw RuntimeError("BaryonPerturbation: invalid baryonFraction < 0 or > 1.");
    }
    if(hubbleConstant <= 0) {
        throw RuntimeError("BaryonPerturbation: invalid hubbleConstant < 0.");
    }
    if(cmbTemperature < 2.7 || cmbTemperature > 2.8) {
        throw RuntimeError("BaryonPerturbation: unexpected cmbTemperature < 2.7 or > 2.8.");
    }
}

local::BaryonPerturbations::~BaryonPerturbations() { }
