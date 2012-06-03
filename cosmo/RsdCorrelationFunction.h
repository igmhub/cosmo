// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_RSD_CORRELATION_FUNCTION
#define COSMO_RSD_CORRELATION_FUNCTION

#include "cosmo/types.h"

namespace cosmo {
    // Calculates the distorted correlation function in redshift space (r,mu)
    // based on the specified ell=0,2,4 (mono,quad,hexa) input correlation
    // functions xi(ell,r). Anisotropic distortions are calculated in linear
    // perturbation theory, with the distant observer assumption. The original
    // reference for this is Kaiser, MNRAS 227, 1 (1987). See also
    // http://astro.berkeley.edu/~mwhite/teachdir/mini_red_mjw.pdf
	class RsdCorrelationFunction {
	public:
	    // Creates a distorted correlation function from the inputs provided.
		RsdCorrelationFunction(CorrelationFunctionPtr xi0,
		    CorrelationFunctionPtr xi2, CorrelationFunctionPtr xi4);
		virtual ~RsdCorrelationFunction();
		// Sets the value of the redshift distortion parameters to use. If beta2=0,
		// then beta2=beta1 is assumed.
        void setDistortion(double beta1, double beta2 = 0);
        // Evaluates the distorted correlation function at the specified pair average co-moving
        // line-of-sight separation rMpch (in Mpc/h) and mu, the cosine of the angle
        // between the pair separation vector and the line of sight.
        double operator()(double rMpch, double mu) const;
        // Evaluates the undistorted correlation function at the specified pair average co-moving
        // line-of-sight separation rMpch (in Mpc/h), for the specified multipole.
        double operator()(double rMpch, Multipole multipole) const;
	private:
        CorrelationFunctionPtr _xi0, _xi2, _xi4;
        double _C0, _C2, _C4;
	}; // RsdCorrelationFunction
} // cosmo

#endif // COSMO_RSD_CORRELATION_FUNCTION
