// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_ABS_HOMOGENEOUS_UNIVERSE
#define COSMO_ABS_HOMOGENEOUS_UNIVERSE

namespace cosmo {
    // Describes a homogenenous and isotropic universe.
	class AbsHomogeneousUniverse {
	public:
		AbsHomogeneousUniverse();
		virtual ~AbsHomogeneousUniverse();
		// Returns the present-day curvature defined as 1 - Omega(0).
        virtual double getCurvature() const = 0;
		// Returns the normalized Hubble function value H(z)/H(0) at the specified
		// redshift z >= 0.
        virtual double getHubbleFunction(double z) const = 0;
        // Returns the comoving line of sight distance in Mpc/h to an emitter with
        // the specified redshift z >= 0.
        virtual double getLineOfSightComovingDistance(double z) const = 0;
        // Returns the comoving transverse distance scale in Mpc/h/rad between two
        // emitters at the same specified redshift z >= 0. Multiply this value by
        // the observed separation angle (rad) to obtain a physical distance in Mpc/h.
        virtual double getTransverseComovingScale(double z) const = 0;
        // Returns the angular diameter distance of an emitter with the specified
        // redshift z >= 0.
        double getAngularDiameterDistance(double z) const;
        // Returns the luminosity distance of an emitter with the specified redshift z >= 0.
        double getLuminosityDistance(double z) const;
        // Returns the lookback time for an emitter at the specified redshift, defined as
        // the difference between the ages of the universe now and when a photon at
        // cosmological redshift z was emitted. Units are secs/h.
        virtual double getLookbackTime(double z) const = 0;
        // Returns the growth function D1(z) for small-scale perturbations in the absence
        // of neutrino free streaming.
        virtual double getGrowthFunction(double z) const = 0;
	private:
	}; // AbsHomogeneousUniverse
	
	inline double AbsHomogeneousUniverse::getAngularDiameterDistance(double z) const {
        return getTransverseComovingScale(z)/(1+z);
	}
	
	inline double AbsHomogeneousUniverse::getLuminosityDistance(double z) const {
        return getTransverseComovingScale(z)*(1+z);
	}
	
	// Returns the present-day Hubble length c/H0 in Mpc/h
	inline double hubbleLength() {
        return 299792458e-5;
	}
	
	// Returns the present-day Hubble time 1/H0 in secs/h
	inline double hubbleTime() {
        return 3.08568025e17;
	}
	
} // cosmo

#endif // COSMO_ABS_HOMOGENEOUS_UNIVERSE
