// Created 09-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "cosmo/RuntimeError.h"
#include "cosmo/AbsHomogeneousUniverse.h"

#include "likely/Integrator.h"

#include "boost/bind.hpp"
#include "boost/ref.hpp"

#include <cmath>

namespace local = cosmo;

local::TransferFunctionPowerSpectrum::TransferFunctionPowerSpectrum(
TransferFunctionPtr transferFunction, double spectralIndex, double deltaH)
: _transferFunction(transferFunction)
{
    setSpectralIndex(spectralIndex);
    setDeltaH(deltaH);
}

local::TransferFunctionPowerSpectrum::~TransferFunctionPowerSpectrum() { }


void local::TransferFunctionPowerSpectrum::setSpectralIndex(double value) {
    _spectralIndex = value;
}

void local::TransferFunctionPowerSpectrum::setDeltaH(double value) {
    if(value <= 0) {
        throw RuntimeError("TransferFunctionPowerSpectrum: invalid deltaH <= 0.");
    }
    _deltaH = value;
    _deltaHSq = value*value;
}

double local::TransferFunctionPowerSpectrum::operator()(double kMpch) const {
    double Tf((*_transferFunction)(kMpch));
    double y(kMpch*hubbleLength());
    return _deltaHSq*std::pow(y,3+_spectralIndex)*Tf*Tf;
}

double local::TransferFunctionPowerSpectrum::setSigma(double sigma, double rMpch, bool gaussian) {
    PowerSpectrumPtr self(new cosmo::PowerSpectrum(boost::ref(*this)));
    double sigmaOld = getRmsAmplitude(self, rMpch, gaussian);
    double ratio(sigma/sigmaOld);
    _deltaH *= ratio;
    _deltaHSq = _deltaH*_deltaH;
    return sigmaOld;
}

namespace cosmo {
    class RmsIntegrand {
    public:
        RmsIntegrand(PowerSpectrumPtr powerSpectrum, double rMpch, bool gaussian)
        : _powerSpectrum(powerSpectrum), _rMpch(rMpch), _gaussian(gaussian) { }
        double operator()(double kMpch) const {
            double kr(kMpch*_rMpch), kr2(kr*kr);
            double wgt(_gaussian ?
                std::exp(-kr2/2) : (std::sin(kr)-kr*std::cos(kr))*3/(kr2*kr));
            return (*_powerSpectrum)(kMpch)*wgt*wgt/kMpch;
        }
    private:
        PowerSpectrumPtr _powerSpectrum;
        double _rMpch;
        bool _gaussian;
    };
} // cosmo::

double local::getRmsAmplitude(PowerSpectrumPtr powerSpectrum, double rMpch, bool gaussian) {
    RmsIntegrand rmsIntegrand(powerSpectrum,rMpch,gaussian);
    likely::Integrator::IntegrandPtr integrand(new likely::Integrator::Integrand(
        boost::ref(rmsIntegrand)));
    likely::Integrator integrator(integrand,1e-8,1e-6);
    return std::sqrt(integrator.integrateSingular(0,1) + integrator.integrateUp(1));
}

template <class P>
local::GenericFunctionPtr local::createFunctionPtr(boost::shared_ptr<P> pimpl) {
    GenericFunctionPtr fptr(new GenericFunction(boost::bind(&P::operator(),pimpl,_1)));
    return fptr;
}

// explicit template instantiations

#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "cosmo/BroadbandPower.h"
#include "likely/Interpolator.h"

template local::PowerSpectrumPtr local::createFunctionPtr<local::TransferFunctionPowerSpectrum>
    (boost::shared_ptr<TransferFunctionPowerSpectrum> pimpl);
template local::PowerSpectrumPtr local::createFunctionPtr<local::BroadbandPower>
    (boost::shared_ptr<BroadbandPower> pimpl);
template local::GenericFunctionPtr local::createFunctionPtr<likely::Interpolator>
    (boost::shared_ptr<likely::Interpolator> pimpl);
