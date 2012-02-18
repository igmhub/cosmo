// Created 18-Feb-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/BroadbandPower.h"
#include "cosmo/RuntimeError.h"

namespace local = cosmo;

local::BroadbandPower::BroadbandPower(int n0, std::vector<double> coefs,
double rmin, double rmax, double r0, double sigmaSq)
: _n0(n0), _coefs(coefs), _rmin(rmin), _rmax(rmax)
{
    if(rmax <= rmin) throw RuntimeError("BroadbandPower: expected rmax > rmin.");
    if(r0 > 0) {
        if(sigmaSq <= 0) throw RuntimeError("BroadbandPower: expected sigmaSq > 0.");
    }
}

local::BroadbandPower::~BroadbandPower() { }
