// Created 09-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/PowerSpectrumCorrelationFunction.h"
#include "cosmo/RuntimeError.h"

namespace local = cosmo;

local::PowerSpectrumCorrelationFunction::PowerSpectrumCorrelationFunction(
PowerSpectrumPtr powerSpectrum, double rmin, double rmax, int nr)
: _powerSpectrum(powerSpectrum), _rmin(rmin), _rmax(rmax), _nr(nr)
{
    if(rmin < 0) {
        throw RuntimeError("PowerSpectrumCorrelationFunction: invalid rmin < 0.");
    }
    if(rmax <= rmin) {
        throw RuntimeError("PowerSpectrumCorrelationFunction: invalid rmax <= rmin.");
    }
    if(nr < 2) {
        throw RuntimeError("PowerSpectrumCorrelationFunction: invalid nr < 2.");
    }
}

local::PowerSpectrumCorrelationFunction::~PowerSpectrumCorrelationFunction() { }

double local::PowerSpectrumCorrelationFunction::operator()(double rMpch) const {
}