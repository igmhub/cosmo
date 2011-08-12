// Created 12-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/FftGaussianRandomFieldGenerator.h"

namespace local = cosmo;

local::FftGaussianRandomFieldGenerator::FftGaussianRandomFieldGenerator(
PowerSpectrumPtr powerSpectrum, double spacing, int nx, int ny, int nz)
: AbsGaussianRandomFieldGenerator(powerSpectrum,spacing,nx,ny,nz)
{
}

local::FftGaussianRandomFieldGenerator::~FftGaussianRandomFieldGenerator() { }

void local::FftGaussianRandomFieldGenerator::_generate(int seed) {
    
}

double local::FftGaussianRandomFieldGenerator::_getFieldUnchecked(int x, int y, int z) const {
    return 0;
}
