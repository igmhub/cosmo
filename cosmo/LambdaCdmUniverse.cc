// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/LambdaCdmUniverse.h"
#include "cosmo/RuntimeError.h"

#include <cmath>

namespace local = cosmo;

local::LambdaCdmUniverse::LambdaCdmUniverse(double OmegaLambda, double OmegaMatter,
double zmax, int nz, double epsAbs)
: _OmegaLambda(OmegaLambda), _OmegaMatter(OmegaMatter), _curvature(1-OmegaLambda-OmegaMatter),
HomogeneousUniverseCalculator(zmax,nz,epsAbs)
{
    if(OmegaLambda < 0) {
        throw RuntimeError("LambdaCdmUniverse: invalid OmegaLambda < 0.");
    }
    if(OmegaMatter < 0) {
        throw RuntimeError("LambdaCdmUniverse: invalid OmegaMatter < 0.");
    }
}

local::LambdaCdmUniverse::~LambdaCdmUniverse() { }

double local::LambdaCdmUniverse::getCurvature() const {
    return _curvature;
}

double local::LambdaCdmUniverse::getHubbleFunction(double z) const {
    if(z < 0) {
        throw RuntimeError("LambdaCdmUniverse::getHubbleFunction: z < 0.");
    }
    // Argument of sqrt is always >= 0 for z >= 0 and each Omega >= 0.
    double ainv(1+z);
    return std::sqrt(_OmegaLambda + ainv*ainv*(_curvature + ainv*_OmegaMatter));
}
