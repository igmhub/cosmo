// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/HomogeneousUniverseCalculator.h"
#include "cosmo/RuntimeError.h"

#include "likely/Integrator.h"
#include "likely/Interpolator.h"

#include "boost/bind.hpp"

#include <cmath>

namespace local = cosmo;

local::HomogeneousUniverseCalculator::HomogeneousUniverseCalculator(double zmax, int nz,
double epsAbs)
: _zmax(zmax), _nz(nz), _epsAbs(epsAbs), _curvatureScale(0)
{
    if(zmax <= 0) {
        throw RuntimeError("HomogeneousUniverseCalculator: invalid zmax <= 0.");
    }
    if(nz < 2) {
        throw RuntimeError("HomogeneousUniverseCalculator: invalid nz < 2.");        
    }
}

local::HomogeneousUniverseCalculator::~HomogeneousUniverseCalculator() { }

double local::HomogeneousUniverseCalculator::getLineOfSightComovingDistance(double z) const {    
    if(z > _zmax) {
        throw RuntimeError("BasicCosmology::getLineOfSightComovingDistance: z > zmax.");
    }
    if(z < 0) {
        throw RuntimeError("BasicCosmology::getLineOfSightComovingDistance: z < 0.");        
    }
    if(!_lineOfSightInterpolator) {
        // Create the interpolator the first time we are called.
        likely::Integrator::IntegrandPtr integrand(new likely::Integrator::Integrand(
            boost::bind(&HomogeneousUniverseCalculator::_lineOfSightIntegrand,this,_1)));
        likely::Integrator integrator(integrand,_epsAbs,0);
        likely::Interpolator::CoordinateValues zValues(_nz), fValues(_nz);
        double dz = _zmax/(_nz-1);
        zValues[0] = 0;
        fValues[0] = 0;
        zValues[1] = dz;
        fValues[1] = integrator.integrateSmooth(0,dz);
        for(int i = 2; i < _nz; ++i) {
            zValues[i] = i*dz;
            fValues[i] = fValues[i-1] + integrator.integrateSmooth((i-1)*dz,i*dz);
        }
        _lineOfSightInterpolator.reset(new likely::Interpolator(zValues,fValues,"cspline"));
    }
    return (*_lineOfSightInterpolator)(z);
}

double local::HomogeneousUniverseCalculator::_lineOfSightIntegrand(double z) const {
    return hubbleLength()/getHubbleFunction(z);
}

double local::HomogeneousUniverseCalculator::getTransverseComovingScale(double z) const {
    double OmegaK(getCurvature());
    double lineOfSight(getLineOfSightComovingDistance(z));
    // All done if there is no curvature.
    if(0 == OmegaK) return lineOfSight;
    // Calculate and remember the value of (c/H0)/sqrt(|curvature|) if necessary.
    if(0 == _curvatureScale) {
        _curvatureScale = hubbleLength()/std::sqrt(std::fabs(OmegaK));
    }
    // Correct the line of sight scale for the curvature.
    if(OmegaK > 0) {
        return _curvatureScale*std::sinh(lineOfSight/_curvatureScale);
    }
    else {
        return _curvatureScale*std::sin(lineOfSight/_curvatureScale);
    }
}

double local::HomogeneousUniverseCalculator::getGrowthFunction(double z) const {
    if(z > _zmax) {
        throw RuntimeError("BasicCosmology::getGrowthFunction: z > zmax.");
    }
    if(z < 0) {
        throw RuntimeError("BasicCosmology::getGrowthFunction: z < 0.");        
    }
    if(!_growthInterpolator) {
        // Create the interpolator the first time we are called.
        likely::Integrator::IntegrandPtr integrand(new likely::Integrator::Integrand(
            boost::bind(&HomogeneousUniverseCalculator::_growthIntegrand,this,_1)));
        likely::Integrator integrator(integrand,_epsAbs,0);
        likely::Interpolator::CoordinateValues zValues(_nz), fValues(_nz);
        double dz = _zmax/(_nz-1);
        zValues[_nz-1] = _zmax;
        fValues[_nz-1] = integrator.integrateUp(_zmax);
        for(int i = _nz-2; i >= 0; --i) {
            double z(i*dz);
            zValues[i] = z;
            fValues[i] = fValues[i+1] + integrator.integrateSmooth(z,(i+1)*dz);
        }
        for(int i = 0; i < _nz; ++i) {
            fValues[i] *= getHubbleFunction(i*dz);
        }
        _growthInterpolator.reset(new likely::Interpolator(zValues,fValues,"cspline"));
    }
    return (*_growthInterpolator)(z);
}

double local::HomogeneousUniverseCalculator::_growthIntegrand(double z) const {
    double hz(getHubbleFunction(z));
    return (1+z)/(hz*hz*hz);
}

double local::HomogeneousUniverseCalculator::getLookbackTime(double z) const {
    if(z > _zmax) {
        throw RuntimeError("BasicCosmology::getLookbackTime: z > zmax.");
    }
    if(z < 0) {
        throw RuntimeError("BasicCosmology::getLookbackTime: z < 0.");        
    }
    if(!_lookbackInterpolator) {
        // Create the interpolator the first time we are called.
        likely::Integrator::IntegrandPtr integrand(new likely::Integrator::Integrand(
            boost::bind(&HomogeneousUniverseCalculator::_lookbackIntegrand,this,_1)));
        likely::Integrator integrator(integrand,_epsAbs,1e-8); // needs epsRel > 0
        likely::Interpolator::CoordinateValues zValues(_nz), fValues(_nz);
        double dz = _zmax/(_nz-1);
        zValues[0] = 0;
        fValues[0] = 0;
        zValues[1] = dz;
        fValues[1] = integrator.integrateSmooth(0,dz);
        for(int i = 2; i < _nz; ++i) {
            zValues[i] = i*dz;
            fValues[i] = fValues[i-1] + integrator.integrateSmooth((i-1)*dz,i*dz);
        }
        _lookbackInterpolator.reset(new likely::Interpolator(zValues,fValues,"cspline"));
    }
    return (*_lookbackInterpolator)(z);    
}

double local::HomogeneousUniverseCalculator::_lookbackIntegrand(double z) const {
    return hubbleTime()/(1+z)/getHubbleFunction(z);
}
