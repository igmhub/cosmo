// Created 18-Feb-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_BROADBAND_POWER
#define COSMO_BROADBAND_POWER

namespace cosmo {
    // Represents a broadband power spectrum P(k) proportional to (1/k)^p.
	class BroadbandPower {
	public:
	    // Creates a new broadband model P(k) = coef*PB(k,p) where PB(k,p)
	    // equals B(p)*(1/k)^p for kmin << k << 1/rmin, with additional terms to
	    // regulate 3D Fourier integrals beyond this range:
	    //
	    //   PB(k,p) = B(p) exp[-(k rmin)^2] / (kmin^p + k^p)
	    //
	    // Parameters kmin and rmin should be specified in h/Mpc and Mpc/h,
	    // respectively, and either or both can be zero to disable the corresponding
	    // regulator term above.
	    // If r0 equals zero (the default), then B(p) = 1. A value r0 > 0 specifies
	    // a scale in Mpc/h used to fix B(p) so that the fluctuations of PB(k,p)
	    // within a top-hat window of radius r0 have a variance sigmaSq. This can
	    // be useful to set a "natural" relative normalization for each term but,
	    // in general, the values B(p) will depend on the choice of kmin, rmin.
		BroadbandPower(double coef, double p,
		    double kmin, double rmin, double r0 = 0, double sigmaSq = 0);
		virtual ~BroadbandPower();
        // Returns the value of k^3/(2pi^2) P(k) for an input wavenumber k in h/Mpc.
        double operator()(double kMpch) const;
	private:
        double _coef, _p, _kmin, _rmin, _kminp, _twopi2;
	}; // BroadbandPower
} // cosmo

#endif // COSMO_BROADBAND_POWER
