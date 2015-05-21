// Created 11-May-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>
// A driver and test program for the DistortedPowerCorrelationHybrid class.
// Calculates the 3D correlation function corresponding to a distorted power spectrum.

#include "cosmo/cosmo.h"
#include "likely/likely.h"
#include "likely/function_impl.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/bind.hpp"
#include "boost/lexical_cast.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

namespace po = boost::program_options;
namespace lk = likely;

class LyaDistortion {
// A simple distortion model for autocorrelations, including linear redshift space effects
// (bias,beta), non-linear large-scale broadening (snlPar,snlPerp), radiation effects
// (biasGamma,biasSourceAbsorber,biasAbsorberResponse,meanFreePath), continuum fitting
// broadband distortion (kc,pc), and non-linear correction (qnl,kv,av,bv,kp)
// or alternative non-linear correction (knl,pnl,kpp,pp,kv0,pv,kvi,pvi).
public:
    LyaDistortion(double bias, double biasbeta,
        double biasGamma, double biasSourceAbsorber, double biasAbsorberResponse,
        double meanFreePath, double snlPar, double snlPerp, double kc, double kcAlt,
        double pc, double sigma8, double qnl, double kv, double av, double bv, double kp,
        double knl, double pnl, double kpp, double pp, double kv0, double pv,
        double kvi, double pvi) :
        _bias(bias), _biasbeta(biasbeta),
        _biasGamma(biasGamma), _biasSourceAbsorber(biasSourceAbsorber),
        _biasAbsorberResponse(biasAbsorberResponse), _meanFreePath(meanFreePath),
        _snlPar2(snlPar*snlPar), _snlPerp2(snlPerp*snlPerp), _kc(kc), _kcAlt(kcAlt), _pc(pc),
        _sigma8(sigma8), _qnl(qnl), _kv(kv), _av(av), _bv(bv), _kp(kp), _knl(knl), _pnl(pnl),
        _kpp(kpp), _pp(pp), _kv0(kv0), _pv(pv), _kvi(kvi), _pvi(pvi)
    {
        _radStrength = biasGamma*biasSourceAbsorber;
    }
    double operator()(double k, double mu, double pk) const {
        // Calculate the k-dependent effective bias
        double bias(_bias);
        if(_radStrength != 0) {
            double s(_meanFreePath*k);
            double Ws = std::atan(s)/s;
            bias += _radStrength*Ws/(1+_biasAbsorberResponse*Ws);
        }
        // Calculate the k-dependent effective beta
        double beta(_biasbeta/bias);
        // Calculate the overall large-scale Lya tracer bias
        double mu2(mu*mu);
        double linear = bias*(1 + beta*mu2);
        // Calculate non-linear broadening
        double snl2 = _snlPar2*mu2 + (1 - mu2)*_snlPerp2;
        double nonlinear = std::exp(-0.5*k*k*snl2);
        // Calculate continuum fitting distortion
        double kpar = std::fabs(k*mu);
        double contdistortion(1);
        if(_kc != 0) {
        	double k1 = std::pow(kpar/_kc + 1,0.75);
    		contdistortion = std::pow((k1-1/k1)/(k1+1/k1),_pc);
        }
        if(_kcAlt != 0) {
        	contdistortion = std::tanh(std::pow(kpar/_kcAlt,_pc));
        }
        // Calculate non-linear correction
        double growth, pecvelocity, pressure, nlcorrection(1);
        if(_qnl != 0) {
        	double sigma8Sim(0.8338);
        	double pi(4*std::atan(1));
        	pk = pk*(sigma8Sim/_sigma8)*(sigma8Sim/_sigma8);
        	double dk = k*k*k*pk/(2*pi*pi);
        	growth = _qnl*dk;
        	pecvelocity = std::pow(k/_kv,_av)*std::pow(std::fabs(mu),_bv);
        	pressure = (k/_kp)*(k/_kp);
        	nlcorrection = std::exp(growth*(1-pecvelocity)-pressure);
        }
        // Calculate alternative non-linear correction (McDonald 2003)
        if(_knl != 0) {
        	double kvel = _kv0*std::pow(1+k/_kvi,_pvi);
        	growth = std::pow(k/_knl,_pnl);
        	pressure = std::pow(k/_kpp,_pp);
        	pecvelocity = std::pow(kpar/kvel,_pv);
        	nlcorrection = std::exp(growth-pressure-pecvelocity);
        }
        return contdistortion*nonlinear*nlcorrection*linear*linear;
    }
private:
    double _bias,_biasbeta,_biasGamma,_biasSourceAbsorber,_biasAbsorberResponse,
    	_meanFreePath, _snlPar2,_snlPerp2,_kc,_kcAlt,_pc,_sigma8,_qnl,_kv,_av,_bv,_kp,
    	_knl,_pnl,_kpp,_pp,_kv0,_pv,_kvi,_pvi,_radStrength;
};

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology distorted power correlation function");
    std::string input,delta,output;
    int nx,ny,nr,nk,nmu;
    double kxmax,spacing,epsAbs,rmin,rmax,maxRelError,kmin,kmax;
    double bias,biasbeta,biasGamma,biasSourceAbsorber,biasAbsorberResponse,meanFreePath,
        snlPar,snlPerp,kc,kcAlt,pc,sigma8,qnl,kv,av,bv,kp,knl,pnl,kpp,pp,kv0,pv,kvi,pvi;
    cli.add_options()
        ("help,h", "prints this info and exits.")
        ("verbose", "prints additional information.")
        ("input,i", po::value<std::string>(&input)->default_value(""),
            "filename to read k,P(k) values from")
        ("delta", po::value<std::string>(&delta)->default_value(""),
            "optional filename of k,P(k) values to subtract from input")
        ("output,o", po::value<std::string>(&output)->default_value(""),
            "base name for saving results")
        ("kxmax", po::value<double>(&kxmax)->default_value(4),
            "Maximum value along x-axis in h/Mpc.")
        ("nx", po::value<int>(&nx)->default_value(400),
            "Grid size along x-axis.")
        ("spacing", po::value<double>(&spacing)->default_value(4),
            "Grid spacing along line-of-sight y-axis in Mpc/h.")
        ("ny", po::value<int>(&ny)->default_value(400),
            "Grid size along line-of-sight y-axis.")
        ("epsAbs", po::value<double>(&epsAbs)->default_value(1e-6),
            "maximum allowed absolute error for 1D integral")
        ("bias", po::value<double>(&bias)->default_value(-0.14),
            "linear tracer bias")
        ("biasbeta", po::value<double>(&biasbeta)->default_value(-0.196),
            "product of bias and linear redshift-space distortion parameter beta")
        ("bias-gamma", po::value<double>(&biasGamma)->default_value(0),
            "Lyman-alpha photoionization bias (nominally 0.13; zero for no UV fluctuations)")
        ("bias-source-absorber", po::value<double>(&biasSourceAbsorber)->default_value(1.0),
            "Lyman-alpha source - absorber bias difference")
        ("bias-absorber-response", po::value<double>(&biasAbsorberResponse)->default_value(-2./3.),
            "Lyman-alpha absorber response bias")
        ("mean-free-path", po::value<double>(&meanFreePath)->default_value(300.),
            "effective mean free path of an ionizing photon in Mpc/h")
        ("snl-par", po::value<double>(&snlPar)->default_value(0.),
            "parallel component of non-linear broadening in Mpc/h")
        ("snl-perp", po::value<double>(&snlPerp)->default_value(0.),
            "perpendicular component of non-linear broadening in Mpc/h")
        ("kc", po::value<double>(&kc)->default_value(0.),
            "parameter for broadband distortion model in h/Mpc (or zero for no distortion)")
        ("kcAlt", po::value<double>(&kcAlt)->default_value(0.),
            "parameter for alternative broadband distortion model in h/Mpc (or zero for no distortion)")
        ("pc", po::value<double>(&pc)->default_value(1.),
            "exponent for broadband distortion model (ignored when kc = 0 and kcAlt = 0)")
        ("sigma8", po::value<double>(&sigma8)->default_value(0.8338),
            "Amplitude of the linear matter power spectrum on the scale of 8 Mpc/h")
        ("qnl", po::value<double>(&qnl)->default_value(0.),
            "strength of non-linear growth correction (nominally 0.867; zero for no non-linear correction)")
        ("kv", po::value<double>(&kv)->default_value(1.05),
            "scale for non-linear peculiar velocity correction in h/Mpc (ignored when qnl = 0)")
        ("av", po::value<double>(&av)->default_value(0.514),
            "exponent for non-linear peculiar velocity correction (ignored when qnl = 0)")
        ("bv", po::value<double>(&bv)->default_value(1.60),
            "exponent for angular dependence of non-linear peculiar velocity correction (ignored when qnl = 0)")
        ("kp", po::value<double>(&kp)->default_value(19.3),
            "scale for non-linear pressure correction in h/Mpc (ignored when qnl = 0)")
        ("knl", po::value<double>(&knl)->default_value(0.),
            "scale for alternative non-linear growth correction in h/Mpc (nominally 6.4; zero for no non-linear correction)")
        ("pnl", po::value<double>(&pnl)->default_value(0.569),
            "exponent for alternative non-linear growth correction (ignored when knl = 0)")
        ("kpp", po::value<double>(&kpp)->default_value(15.3),
            "scale for alternative non-linear pressure correction in h/Mpc (ignored when knl = 0)")
        ("pp", po::value<double>(&pp)->default_value(2.01),
            "exponent for alternative non-linear pressure correction (ignored when knl = 0)")
        ("kv0", po::value<double>(&kv0)->default_value(1.22),
            "scale for alternative line-of-sight non-linear peculiar velocity correction in h/Mpc (ignored when knl = 0)")
        ("pv", po::value<double>(&pv)->default_value(1.5),
            "exponent for alternative line-of-sight non-linear peculiar velocity correction (ignored when knl = 0)")
        ("kvi", po::value<double>(&kvi)->default_value(0.923),
            "scale for alternative isotropic non-linear peculiar velocity correction in h/Mpc (ignored when knl = 0)")
        ("pvi", po::value<double>(&pvi)->default_value(0.451),
            "exponent for alternative isotropic non-linear peculiar velocity correction (ignored when knl = 0)")
        ("max-rel-error", po::value<double>(&maxRelError)->default_value(1e-3),
            "maximum allowed relative error for power-law extrapolation of input P(k)")
        ("kmin", po::value<double>(&kmin)->default_value(0.0001),
            "minimum value of comoving separation to use")
        ("kmax", po::value<double>(&kmax)->default_value(1152.5),
            "maximum value of comoving separation to use")
        ("nk", po::value<int>(&nk)->default_value(814),
            "number of log-spaced k values for saving results")
        ("rmin", po::value<double>(&rmin)->default_value(10.),
            "minimum value of comoving separation to use")
        ("rmax", po::value<double>(&rmax)->default_value(200.),
            "maximum value of comoving separation to use")
        ("nr", po::value<int>(&nr)->default_value(191),
            "number of points spanning [rmin,rmax] to use")
        ("nmu", po::value<int>(&nmu)->default_value(11),
            "number of equally spaced mu_k and mu_r values for saving results")
        ;
    // Do the command line parsing now
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose"));

    if(input.length() == 0) {
        std::cerr << "Missing input filename." << std::endl;
        return 1;
    }
    // Calculate maximum r value to use for bicubic interpolation grid.
    double rgridmax = rmax + spacing;

    try {
        cosmo::TabulatedPowerCPtr power =
            cosmo::createTabulatedPower(input,true,true,maxRelError,verbose);
        if(delta.length() > 0) {
            cosmo::TabulatedPowerCPtr power2 =
                cosmo::createTabulatedPower(delta,true,true,maxRelError,verbose);
            power = power->createDelta(power2,verbose);
        }
        lk::GenericFunctionPtr PkPtr =
            lk::createFunctionPtr<const cosmo::TabulatedPower>(power);

        boost::shared_ptr<LyaDistortion> rsd(new LyaDistortion(
            bias,biasbeta,biasGamma,biasSourceAbsorber,biasAbsorberResponse,meanFreePath,
            snlPar,snlPerp,kc,kcAlt,pc,sigma8,qnl,kv,av,bv,kp,knl,pnl,kpp,pp,kv0,pv,kvi,pvi));
        cosmo::KMuPkFunctionCPtr distPtr(new cosmo::KMuPkFunction(boost::bind(
            &LyaDistortion::operator(),rsd,_1,_2,_3)));

        double kxmin = power->getKMin();
    	cosmo::DistortedPowerCorrelationHybrid dpc(PkPtr,distPtr,kxmin,kxmax,nx,spacing,ny,rgridmax,epsAbs);
    	if(verbose) {
        	std::cout << "Memory size = "
            	<< boost::format("%.1f Mb") % (dpc.getMemorySize()/1048576.) << std::endl;
    	}
    	// Transform
    	dpc.transform();
        if(output.length() > 0) {
            double dmu = 1./(nmu-1.);
            // Write out values tabulated for log-spaced k
            double dk = std::pow(kmax/kmin,1./(nk-1.));
            std::string kfile = output + ".k.dat";
            std::ofstream kout(kfile.c_str());
            for(int i = 0; i < nk; ++i) {
                double k = kmin*std::pow(dk,i);
                kout << boost::lexical_cast<std::string>(k);
                for(int j = 0; j < nmu; ++j) {
                    double mu = 1 - j*dmu;
                    kout << ' ' << boost::lexical_cast<std::string>(dpc.getPower(k,mu));
                }
                kout << std::endl;
            }
            kout.close();
            // Write out values tabulated for linear-spaced r
            double dr = (rmax-rmin)/(nr-1.);
            std::string rfile = output + ".r.dat";
            std::ofstream rout(rfile.c_str());
            for(int i = 0; i < nr; ++i) {
                double r = rmin + i*dr;
                rout << boost::lexical_cast<std::string>(r);
                for(int j = 0; j < nmu; ++j) {
                    double mu = 1 - j*dmu;
                    rout << ' ' << boost::lexical_cast<std::string>(dpc.getCorrelation(r,mu));
                }
                rout << std::endl;                
            }
            rout.close();
        }
    }
    catch(std::runtime_error const &e) {
        std::cerr << "ERROR: exiting with an exception:\n  " << e.what() << std::endl;
        return -1;
    }

    return 0;
}