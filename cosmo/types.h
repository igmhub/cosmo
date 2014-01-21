// Created 8-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_TYPES
#define COSMO_TYPES

#include "boost/function.hpp"
#include "boost/smart_ptr.hpp"

namespace cosmo {
    
    enum Multipole { Monopole = 0, Quadrupole = 2, Hexadecapole = 4 };
    
    class AbsHomogeneousUniverse;
    typedef boost::shared_ptr<AbsHomogeneousUniverse> AbsHomogeneousUniversePtr;
    
    class RsdCorrelationFunction;
    typedef boost::shared_ptr<RsdCorrelationFunction> RsdCorrelationFunctionPtr;

    class AbsGaussianRandomFieldGenerator;
    typedef boost::shared_ptr<AbsGaussianRandomFieldGenerator> AbsGaussianRandomFieldGeneratorPtr;

    class TabulatedPower;
    typedef boost::shared_ptr<const TabulatedPower> TabulatedPowerCPtr;

    // Represents a function that returns a dimensionless transfer function value T(k)
    // given an input wavenumber k in 1/(Mpc/h).
    typedef boost::function<double (double)> TransferFunction;
    typedef boost::shared_ptr<TransferFunction> TransferFunctionPtr;
    
    // Represents a function that returns the power spectrum value k^3/(2pi^2) P(k)
    // given an input wavenumber k in 1/(Mpc/h).
    typedef boost::function<double (double)> PowerSpectrum;
    typedef boost::shared_ptr<PowerSpectrum> PowerSpectrumPtr;
    
    // Represents a function that returns a correlation function value given an
    // input pair separation in Mpc/h.
    typedef boost::function<double (double)> CorrelationFunction;
    typedef boost::shared_ptr<CorrelationFunction> CorrelationFunctionPtr;

    // Represents a function of (r,mu)
    typedef boost::function<double (double,double)> RMuFunction;
    typedef boost::shared_ptr<const RMuFunction> RMuFunctionCPtr;
    
} // cosmo

#endif // COSMO_TYPES
