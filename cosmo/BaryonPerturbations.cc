// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/BaryonPerturbations.h"
#include "cosmo/RuntimeError.h"
#include "cosmo/AbsHomogeneousUniverse.h"

#include <cmath>

namespace local = cosmo;

local::BaryonPerturbations::BaryonPerturbations(double omegaMatter, double omegaBaryon,
double hubbleConstant, double cmbTemperature, BaoOption baoOption)
: _omegaMatter(omegaMatter), _omegaBaryon(omegaBaryon),
_hubbleConstant(hubbleConstant), _cmbTemperature(cmbTemperature), _baoOption(baoOption)
{
    if(omegaMatter < 0) {
        throw RuntimeError("BaryonPerturbation: invalid omegaMatter < 0.");
    }
    if(omegaBaryon < 0 || omegaBaryon > omegaMatter) {
        throw RuntimeError("BaryonPerturbation: invalid omegaBaryon < 0 or > omegaMatter.");
    }
    if(hubbleConstant <= 0) {
        throw RuntimeError("BaryonPerturbation: invalid hubbleConstant < 0.");
    }
    if(cmbTemperature < 2.7 || cmbTemperature > 2.8) {
        throw RuntimeError("BaryonPerturbation: unexpected cmbTemperature < 2.7 or > 2.8.");
    }

    double f_baryon(_omegaBaryon/_omegaMatter);
    double f2(f_baryon*f_baryon),f3(f2*f_baryon);
    
    double hSq(_hubbleConstant*_hubbleConstant);
    _omhh = _omegaMatter*hSq;    
    _obhh = _omegaBaryon*hSq;

    _theta_cmb = _cmbTemperature/2.7;
    double tcmb2 = _theta_cmb*_theta_cmb;
    double tcmb4 = tcmb2*tcmb2;

    _z_equality = 2.50e4*_omhh/tcmb4;  /* Really 1+z */
    _k_equality = 0.0746*_omhh/tcmb2;

    double z_drag_b1(0.313*std::pow(_omhh,-0.419)*(1+0.607*std::pow(_omhh,0.674)));
    double z_drag_b2(0.238*std::pow(_omhh,0.223));
    _z_drag = 1291*std::pow(_omhh,0.251)/(1+0.659*std::pow(_omhh,0.828))*
		(1+z_drag_b1*std::pow(_obhh,z_drag_b2));
    
    _R_drag = 31.5*_obhh/tcmb4*(1000/(1+_z_drag));
    _R_equality = 31.5*_obhh/tcmb4*(1000/_z_equality);

    _sound_horizon = 2./3./_k_equality*std::sqrt(6./_R_equality)*
        std::log((std::sqrt(1+_R_drag)+
	    std::sqrt(_R_drag+_R_equality))/(1+std::sqrt(_R_equality)));
	    
    _k_silk = 1.6*std::pow(_obhh,0.52)*std::pow(_omhh,0.73)*(1+std::pow(10.4*_omhh,-0.95));

    double alpha_c_a1(std::pow(46.9*_omhh,0.670)*(1+std::pow(32.1*_omhh,-0.532)));
    double alpha_c_a2(std::pow(12.0*_omhh,0.424)*(1+std::pow(45.0*_omhh,-0.582)));
    _alpha_c = std::pow(alpha_c_a1,-f_baryon)*std::pow(alpha_c_a2,-f3);
    
    double beta_c_b1(0.944/(1+std::pow(458*_omhh,-0.708)));
    double beta_c_b2(std::pow(0.395*_omhh, -0.0266));
    _beta_c = 1.0/(1+beta_c_b1*(std::pow(1-f_baryon, beta_c_b2)-1));

    double y(_z_equality/(1+_z_drag));
    double ytmp(std::sqrt(1+y));
    double alpha_b_G(y*(-6.*ytmp+(2.+3.*y)*std::log((ytmp+1)/(ytmp-1))));
    _alpha_b = 2.07*_k_equality*_sound_horizon*std::pow(1+_R_drag,-0.75)*alpha_b_G;
    
    _beta_node = 8.41*std::pow(_omhh, 0.435);
    _beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*std::sqrt(std::pow(17.2*_omhh,2.0)+1);

    _k_peak = 2.5*3.14159*(1+0.217*_omhh)/_sound_horizon;
    _sound_horizon_fit = 44.5*std::log(9.83/_omhh)/std::sqrt(1+10.0*std::pow(_obhh,0.75));

    _alpha_gamma = 1-0.328*std::log(431.0*_omhh)*f_baryon + 0.38*std::log(22.3*_omhh)*f2;
}

local::BaryonPerturbations::~BaryonPerturbations() { }

#include <iostream>

double local::BaryonPerturbations::getNode(int n) const {
    if(n <= 0) throw RuntimeError("BaryonsPerturbations::getNode: invalid n <= 0.");
    double pi(4*atan(1));
    double n2(n*n), n3(n2*n), bpi(_beta_node/pi);
    return (_baoOption == ShiftedOscillation) ?
        std::pow((n3+std::sqrt(n3*(n3+4*bpi*bpi*bpi)))/2,1./3.) : n;
}

double local::BaryonPerturbations::getCdmTransfer(double kMpch) const {
    if(kMpch != _kSave) {
        calculateTransferFunctions(kMpch,_Tf_baryon,_Tf_cdm,_Tf_full,_Tf_nw,_baoOption);
        _kSave = kMpch;
    }
    return _Tf_cdm;
}

double local::BaryonPerturbations::getBaryonTransfer(double kMpch) const {
    if(kMpch != _kSave) {
        calculateTransferFunctions(kMpch,_Tf_baryon,_Tf_cdm,_Tf_full,_Tf_nw,_baoOption);
        _kSave = kMpch;
    }
    return _Tf_baryon;
}

double local::BaryonPerturbations::getMatterTransfer(double kMpch) const {
    if(kMpch != _kSave) {
        calculateTransferFunctions(kMpch,_Tf_baryon,_Tf_cdm,_Tf_full,_Tf_nw,_baoOption);
        _kSave = kMpch;
    }
    return _Tf_full;
}

double local::BaryonPerturbations::getNoWigglesTransfer(double kMpch) const {
    if(kMpch != _kSave) {
        calculateTransferFunctions(kMpch,_Tf_baryon,_Tf_cdm,_Tf_full,_Tf_nw,_baoOption);
        _kSave = kMpch;
    }
    return _Tf_nw;
}

void local::BaryonPerturbations::calculateTransferFunctions(double kMpch,
double &Tf_baryon, double &Tf_cdm, double &Tf_full, double &Tf_nw, BaoOption baoOption) const {

    if(0 == kMpch) {
        Tf_baryon = Tf_cdm = Tf_full = 1;
        return;
    }

    double k(kMpch*_hubbleConstant);
    double q(k/13.41/_k_equality);
    double qSq(q*q);
    double xx(k*_sound_horizon);

    double T_c_ln_beta(std::log(2.718282+1.8*_beta_c*q));
    double T_c_ln_nobeta(std::log(2.718282+1.8*q));
    double T_c_C_alpha(14.2/_alpha_c + 386.0/(1+69.9*std::pow(q,1.08)));
    double T_c_C_noalpha(14.2 + 386.0/(1+69.9*std::pow(q,1.08)));

    double tmp(xx/5.4),tmp2(tmp*tmp);
    double T_c_f(1.0/(1.0+tmp2*tmp2));
    Tf_cdm = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*qSq) +
	    (1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*qSq);
    
    double T_b_T0(T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*qSq));

    tmp = xx/5.2;
    Tf_baryon = T_b_T0/(1+tmp*tmp);
    tmp = _beta_b/xx;
    tmp2 = tmp*tmp;
	Tf_baryon += _alpha_b/(1+tmp*tmp2)*std::exp(-std::pow(k/_k_silk,1.4));

    // Add baryon acoustic oscillations
    if(baoOption == PeriodicOscillation) {
        // Force the oscillations to be exactly periodic with the sound horizon scale.
        Tf_baryon *= std::sin(xx)/(xx);
    }
    else {
        // Calculate the effective sound horizon at this wavenumber
	    tmp = _beta_node/xx;
        tmp2 = tmp*tmp;
        double s_tilde(_sound_horizon*std::pow(1+tmp*tmp2,-1./3.));
        double xx_tilde(k*s_tilde);
    	if(baoOption == NoOscillation) {
    	    // See the first paragraph of Section 3.3 in astro-ph/9709112.
    	    // This function roughly follows the central lobe of sin(x)/x,
    	    // then asymptotes to 1/xx_tilde, which is the envelope of sin(x)/x.
            tmp2 = xx_tilde*xx_tilde;
            Tf_baryon *= std::pow(1+tmp2*tmp2,-0.25);
        }
        else {
            Tf_baryon *= std::sin(xx_tilde)/(xx_tilde);
        }
    }
    
    double f_baryon(_obhh/_omhh);
    Tf_full = f_baryon*Tf_baryon + (1-f_baryon)*Tf_cdm;

	// Calculate the no-wiggles transfer function using eqns (29-31) of EH98.
	// Equations below reproduce exactly the sample code from EH98, but should consider
	// replacing _sound_horizon_fit with _sound_horizon and a full precision value for e.
    tmp = 0.43*k*_sound_horizon_fit;
    tmp2 = tmp*tmp;
    double gamma_eff = _omhh*(_alpha_gamma + (1-_alpha_gamma)/(1+tmp2*tmp2));
    double q_eff(q*_omhh/gamma_eff);
    double L0 = std::log(2.0*2.718282+1.8*q_eff); // using e = 2.718282
    double C0 = 14.2 + 731.0/(1+62.5*q_eff);
    Tf_nw = L0/(L0 + C0*q_eff*q_eff);
}
