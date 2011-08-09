// Created 09-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "cosmo/RuntimeError.h"
#include "cosmo/AbsHomogeneousUniverse.h"

#include "likely/Integrator.h"

#include "boost/ref.hpp"

#include <cmath>

namespace local = cosmo;

local::TransferFunctionPowerSpectrum::TransferFunctionPowerSpectrum(
TransferFunctionPtr transferFunction, double spectralIndex, double deltaH)
: _transferFunction(transferFunction), _spectralIndex(spectralIndex), _deltaHSq(deltaH*deltaH)
{
    if(deltaH <= 0) {
        throw RuntimeError("TransferFunctionPowerSpectrum: invalid deltaH <= 0.");
    }
}

local::TransferFunctionPowerSpectrum::~TransferFunctionPowerSpectrum() { }

double local::TransferFunctionPowerSpectrum::operator()(double kMpch) const {
    double Tf((*_transferFunction)(kMpch));
    double y(kMpch*hubbleLength());
    return _deltaHSq*std::pow(y,3+_spectralIndex)*Tf*Tf;
}

namespace cosmo {
    class RmsIntegrand {
    public:
        RmsIntegrand(PowerSpectrumPtr powerSpectrum, double rMpch)
        : _powerSpectrum(powerSpectrum), _rMpch(rMpch) { }
        double operator()(double kMpch) const {
            double kr(kMpch*_rMpch), kr2(kr*kr);
            double wgt((std::sin(kr)-kr*std::cos(kr))*3/(kr2*kr));
            return (*_powerSpectrum)(kMpch)*wgt*wgt/kMpch;
        }
    private:
        PowerSpectrumPtr _powerSpectrum;
        double _rMpch;
    };
} // cosmo::

double local::getRmsAmplitude(PowerSpectrumPtr powerSpectrum, double rMpch) {
    RmsIntegrand rmsIntegrand(powerSpectrum,rMpch);
    likely::Integrator::IntegrandPtr integrand(new likely::Integrator::Integrand(
        boost::ref(rmsIntegrand)));
    likely::Integrator integrator(integrand,1e-8,1e-6);
    return integrator.integrateSingular(0,1) + integrator.integrateUp(1);
}