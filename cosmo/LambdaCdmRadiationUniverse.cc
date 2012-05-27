// Created 26-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/LambdaCdmRadiationUniverse.h"
#include "cosmo/RuntimeError.h"

#include <cmath>

namespace local = cosmo;

local::LambdaCdmRadiationUniverse::LambdaCdmRadiationUniverse(double OmegaMatter, double OmegaK,
    double hValue, double Tcmb, double Nnu, double zmax, int nz, double epsAbs)
: _OmegaMatter(OmegaMatter), _OmegaK(OmegaK), HomogeneousUniverseCalculator(zmax,nz,epsAbs)
{
    if(OmegaMatter < 0) {
        throw RuntimeError("LambdaCdmRadiationUniverse: invalid OmegaMatter < 0.");
    }
    // Calculate H0 in /s
    double H0 = 3.2407764868054625e-18*hValue;
    // Calculate the present critical density in J/m^3. The constant is 3c^2/(8piG).
    double rhoc = 1.6073793277185455e26*H0*H0;
    // Calculate the Stefan-Boltzman energy density for the specified photon temperature.
    // The constant is pi^2/15 / (hbarc)^3 kB^4.
    double rho = 7.565768019130482e-16*std::pow(Tcmb,4);
    // Add massless neutrino energy density. Constant is (7/8) (4/11)^(4/3).
    double fnu = 0.22710731766023895*Nnu;
    rho *= (1+fnu);
    // Calculate the present radiation energy density.
    _OmegaRadiation = rho/rhoc;
    // Calculate the present dark energy density.
    _OmegaLambda = 1 - _OmegaMatter - _OmegaRadiation - _OmegaK;
}

local::LambdaCdmRadiationUniverse::~LambdaCdmRadiationUniverse() { }

double local::LambdaCdmRadiationUniverse::getCurvature() const {
    return _OmegaK;
}

double local::LambdaCdmRadiationUniverse::getHubbleFunction(double z) const {
    if(z < 0) {
        throw RuntimeError("LambdaCdmRadiationUniverse::getHubbleFunction: z < 0.");
    }
    // Argument of sqrt is always >= 0 for z >= 0 and each Omega >= 0.
    double ainv(1+z);
    return std::sqrt(_OmegaLambda + ainv*ainv*(_OmegaK + ainv*(_OmegaMatter + ainv*_OmegaRadiation)));
}
