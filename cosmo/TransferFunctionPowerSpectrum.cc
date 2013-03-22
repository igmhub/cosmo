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

double local::legendreP(int ell, double mu) {
    double mu2(mu*mu);
    switch(ell) {
        case 0:
        return 1;
        case 2:
        return (-1+mu2*3)/2;
        case 4:
        return (3-mu2*(30-mu2*35))/8;
        case 6:
        return (-5+mu2*(105-mu2*(315-mu2*231)))/16;
        case 8:
        return (35-mu2*(1260-mu2*(6930-mu2*(12012-mu2*6435))))/128;
        case 10:
        return (-63+mu2*(3465-mu2*(30030-mu2*(90090-mu2*(109395-mu2*46189)))))/256;
        case 12:
        return (231-mu2*(18018-mu2*(225225-mu2*(1021020-mu2*(2078505-mu2*(1939938-mu2*676039))))))/1024;
        default:
        return 0;
    }
}

namespace cosmo {
    class MultipoleIntegrand {
    public:
        MultipoleIntegrand(GenericFunctionPtr fOfMuPtr, int ell)
        : _fOfMuPtr(fOfMuPtr), _ell(ell) {
            if(ell < 0 || ell % 2 == 1) throw RuntimeError("MultipoleIntegrand: invalid ell.");
            if(ell > 12) throw RuntimeError("MultipoleIntegrand: ell > 12 not implemented yet.");
            _norm = 0.5*(2*ell+1);
        }
        double operator()(double mu) const {
            double fval = (*_fOfMuPtr)(mu);
            return _norm*legendreP(_ell,mu)*fval;
        }
    private:
        GenericFunctionPtr _fOfMuPtr;
        int _ell;
        double _norm;
    };
} // cosmo::

double local::getMultipole(GenericFunctionPtr fOfMuPtr, int ell, double epsAbs, double epsRel) {
    MultipoleIntegrand multipoleIntegrand(fOfMuPtr,ell);
    likely::Integrator::IntegrandPtr integrand(new likely::Integrator::Integrand(
        boost::ref(multipoleIntegrand)));
    likely::Integrator integrator(integrand,epsAbs,epsRel);
    // Factor of two is because we integrate over 0 < mu < 1 instead of -1 < mu < +1.
    return 2*integrator.integrateSmooth(0,1);
}

// explicit template instantiation for creating a function pointer to a TransferFunctionPowerSpectrum.

#include "likely/function_impl.h"

template local::PowerSpectrumPtr likely::createFunctionPtr<local::TransferFunctionPowerSpectrum>
    (boost::shared_ptr<local::TransferFunctionPowerSpectrum> pimpl);
