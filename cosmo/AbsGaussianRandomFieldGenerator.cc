// Created 12-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/AbsGaussianRandomFieldGenerator.h"
#include "cosmo/RuntimeError.h"

#include "likely/Random.h"

namespace local = cosmo;

local::AbsGaussianRandomFieldGenerator::AbsGaussianRandomFieldGenerator(
PowerSpectrumPtr powerSpectrum, double spacing, int nx, int ny, int nz, likely::RandomPtr random)
: _powerSpectrum(powerSpectrum), _spacing(spacing), _nx(nx), _ny(ny), _nz(nz), _random(random)
{
    if(spacing <= 0) {
        throw RuntimeError("AbsGaussianRandomFieldGenerator: invalid spacing <= 0.");
    }
    if(nx <= 0) {
        throw RuntimeError("AbsGaussianRandomFieldGenerator: invalid nx <= 0.");        
    }
    if(ny <= 0) {
        throw RuntimeError("AbsGaussianRandomFieldGenerator: invalid ny <= 0.");        
    }
    if(nz <= 0) {
        throw RuntimeError("AbsGaussianRandomFieldGenerator: invalid nz <= 0.");        
    }
    if(!_random) _random = likely::Random::instance();
}

local::AbsGaussianRandomFieldGenerator::~AbsGaussianRandomFieldGenerator() { }

double local::AbsGaussianRandomFieldGenerator::getField(int x, int y, int z) const {
    if(x < 0 || x >= _nx) {
        throw RuntimeError("AbsGaussianRandomFieldGenerator: invalid x < 0 or >= nx.");
    }
    if(y < 0 || y >= _ny) {
        throw RuntimeError("AbsGaussianRandomFieldGenerator: invalid y < 0 or >= ny.");
    }
    if(z < 0 || z >= _nz) {
        throw RuntimeError("AbsGaussianRandomFieldGenerator: invalid z < 0 or >= nz.");
    }
    return _getFieldUnchecked(x,y,z);
}

std::size_t local::AbsGaussianRandomFieldGenerator::getMemorySize() const {
    return 0;
}

double local::AbsGaussianRandomFieldGenerator::getPower(double k) const {
    return (*_powerSpectrum)(k);
}
