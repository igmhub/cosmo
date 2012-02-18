// Created 18-Feb-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_BROADBAND_POWER
#define COSMO_BROADBAND_POWER

#include <vector>

namespace cosmo {
    // Represents a broadband power spectrum P(k) modeled as a polynomial in (1/k).
	class BroadbandPower {
	public:
	    // Creates a new broadband model P(k) = Sum[coefs[n]*PB(k,n),{n,n0,n0+N}] where
	    // N is the number of elements in the input coefficient array and PB(k,n)
	    // equals B(n)*k^(-n) for 1/rmax << k << 1/rmin, with additional terms to
	    // regulate 3D Fourier integrals beyond this range:
	    //
	    //   PB(k,n) = B(n) exp(-(k*rmin)^2) rmax^n/(1 + (k*rmax)^2)
	    //
	    // If r0 equals zero (the default), then B(n) = 1. A value r0 > 0 specifies
	    // a scale in Mpc/h used to fix B(n) so that the fluctuations of PB(k,n)
	    // within a top-hat window of radius r0 have a variance sigmaSq. This can
	    // be useful to set a "natural" relative normalization for each term but,
	    // in general, the values B(n) will depend on the choice of rmax (and also
	    // of rmin unless rmin << r0/100). rmin and rmax are in Mpc/h.
		BroadbandPower(int n0, std::vector<double> coefs,
		    double rmin, double rmax, double r0 = 0, double sigmaSq = 0);
		virtual ~BroadbandPower();
	private:
        int _n0;
        std::vector<double> _coefs;
        double _rmin, _rmax;
	}; // BroadbandPower
} // cosmo

#endif // COSMO_BROADBAND_POWER
