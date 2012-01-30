// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_BARYON_PERTURBATIONS
#define COSMO_BARYON_PERTURBATIONS

namespace cosmo {
    // Calculates baryon perturbations to a homogenous universe using the results in
    // Eisenstein & Hu, "Baryonic Features in the Matter Transfer Function", astro-ph/9709112
    // Code adapted from http://background.uchicago.edu/~whu/transfer/transferpage.html
	class BaryonPerturbations {
	public:
        // Options for including baryon acoustic oscillations in the transfer function calculation.
        // See Section 3.2 of astro-ph/9709112 for details.
        enum BaoOption {
            NoOscillation,          // sinc(k*s) replaced by 
            PeriodicOscillation,    // nodes are periodically spaced in k
            ShiftedOscillation      // first few nodes (k*s < 10) are shifted to higher k
        };
		BaryonPerturbations(double omegaMatter, double omegaBaryon,
		    double hubbleConstant, double cmbTemperature, BaoOption baoOption = ShiftedOscillation);
		virtual ~BaryonPerturbations();
		// Returns the redshift of matter-radiation equality. See eqn. (2).
        double getMatterRadiationEqualityRedshift() const;
        // Returns the wavenumber of the particle horizon at matter-radiation equality in
        // 1/(Mpc/h). See eqn. (3).
        double getMatterRadiationEqualityScale() const;
        // Returns the "drag epoch", defined as the redshift at which the baryons are
        // released from the Compton drag of the photons. See eqn. (4).
        double getDragEpoch() const;
        // Returns the sound horizon at the drag epoch in Mpc/h as the comoving distance
        // a wave can travel prior to the drag epoch redshift. See eqn. (6).
        double getSoundHorizon() const;
        // Returns the wavenumber in 1/(Mpc/h) characterizing baryon-photon diffusion.
        double getSilkDampingScale() const;
        // Returns the CDM transfer function value at the specified wavenumber in 1/(Mpc/h).
        double getCdmTransfer(double kMpch) const;
        // Returns the baryon transfer function value at the specified wavenumber
        // in 1/(Mpc/h).
        double getBaryonTransfer(double kMpch) const;
        // Returns the CDM + baryon transfer function value at the specified wavenumber
        // in 1/(Mpc/h).
        double getMatterTransfer(double kMpch) const;
		// Calculates and stores the baryon, CDM, and full (baryon+CDM) transfer functions
		// for the specified input wavenumber k in 1/(Mpc/h).
        void calculateTransferFunctions(double kMpch,
            double &Tf_baryon, double &Tf_cdm, double &Tf_full,
            BaoOption baoOption = ShiftedOscillation) const;
	private:
        double _omegaMatter, _omegaBaryon, _hubbleConstant, _cmbTemperature;
        BaoOption _baoOption;
        double
            _omhh,		    /* Omega_matter*h^2 */
        	_obhh,		    /* Omega_baryon*h^2 */
        	_theta_cmb,	    /* Tcmb in units of 2.7 K */
        	_z_equality,	/* Redshift of matter-radiation equality, really 1+z */
        	_k_equality,	/* Scale of equality, in Mpc^-1 */
        	_z_drag,		/* Redshift of drag epoch */
        	_R_drag,		/* Photon-baryon ratio at drag epoch */
        	_R_equality,	/* Photon-baryon ratio at equality epoch */
        	_sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
        	_k_silk,		/* Silk damping scale, in Mpc^-1 */
        	_alpha_c,	    /* CDM suppression */
        	_beta_c,		/* CDM log shift */
        	_alpha_b,	    /* Baryon suppression */
        	_beta_b,		/* Baryon envelope shift */
        	_beta_node,	    /* Sound horizon shift */
        	_k_peak,		/* Fit to wavenumber of first peak, in Mpc^-1 */
        	_sound_horizon_fit,	/* Fit to sound horizon, in Mpc */
        	_alpha_gamma;	/* Gamma suppression in approximate TF */
        mutable double _Tf_baryon, _Tf_cdm, _Tf_full, _kSave;
	}; // BaryonPerturbations
	
	inline double BaryonPerturbations::getMatterRadiationEqualityRedshift() const {
        return _z_equality - 1;
	}
	inline double BaryonPerturbations::getMatterRadiationEqualityScale() const {
        return _k_equality/_hubbleConstant;
	}
	inline double BaryonPerturbations::getDragEpoch() const {
        return _z_drag;
	}
	inline double BaryonPerturbations::getSoundHorizon() const {
        return _sound_horizon*_hubbleConstant;
	}
	inline double BaryonPerturbations::getSilkDampingScale() const {
        return _k_silk/_hubbleConstant;
	}
	
} // cosmo

#endif // COSMO_BARYON_PERTURBATIONS
