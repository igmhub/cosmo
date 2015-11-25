// Created 22-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A driver and test program for the DistortedPowerCorrelation class.
// Calculates the 3D correlation function corresponding to a distorted power spectrum.

// Sample timing results on OS 10.8 laptop using:
//
// ./cosmodpc -i ../../baofit/models/PlanckWPBestFitLCDM_matterpower.dat OPTIONS
//
//   TIME  OPTIONS
//  127.3s --repeat 10000 --direct-power-multipoles
//  126.5s --repeat 10000 --direct-power-multipoles --optimize
//   83.2s --repeat 10000 --direct-power-multipoles --bypass
//   53.2s --repeat 10000
//   53.4s --repeat 10000 --optimize
//   53.2s --repeat 10000 --bypass

#include "cosmo/cosmo.h"
#include "likely/likely.h"
#include "likely/function_impl.h"

#include "boost/program_options.hpp"
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
// broadband distortion (kc,pc), non-linear correction (knl,pnl,kpp,pp,kv0,pv,kvi,pvi)
// and parallel pixelization smoothing.
public:
    LyaDistortion(double bias, double biasbeta,
        double biasGamma, double biasSourceAbsorber, double biasAbsorberResponse,
        double meanFreePath, double snlPar, double snlPerp, double kc, double kcAlt,
        double pc, double knl, double pnl, double kpp, double pp, double kv0, double pv,
        double kvi, double pvi, double pixPar) :
        _bias(bias), _biasbeta(biasbeta),
        _biasGamma(biasGamma), _biasSourceAbsorber(biasSourceAbsorber),
        _biasAbsorberResponse(biasAbsorberResponse), _meanFreePath(meanFreePath),
        _snlPar2(snlPar*snlPar), _snlPerp2(snlPerp*snlPerp),
        _kc(kc), _kcAlt(kcAlt), _pc(pc), _knl(knl), _pnl(pnl),
        _kpp(kpp), _pp(pp), _kv0(kv0), _pv(pv), _kvi(kvi), _pvi(pvi), _pixPar(pixPar)
    {
        _radStrength = biasGamma*biasSourceAbsorber;
    }
    double operator()(double k, double mu) const {
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
        // Calculate non-linear correction (McDonald 2003)
        double growth, pecvelocity, pressure, nlcorrection(1);
        if(_knl != 0) {
        	double kvel = _kv0*std::pow(1+k/_kvi,_pvi);
        	growth = std::pow(k/_knl,_pnl);
        	pressure = std::pow(k/_kpp,_pp);
        	pecvelocity = std::pow(kpar/kvel,_pv);
        	nlcorrection = std::exp(growth-pressure-pecvelocity);
        }
        // Parallel pixelization smoothing
        double pixelization(1);
        if(_pixPar != 0) {
            double pix = std::sin(_pixPar*kpar)/(_pixPar*kpar);
            pixelization = pix*pix;
        }
        return contdistortion*nonlinear*nlcorrection*linear*linear*pixelization;
    }
private:
    double _bias,_biasbeta,_biasGamma,_biasSourceAbsorber,_biasAbsorberResponse,
        _meanFreePath,_snlPar2,_snlPerp2,_kc,_kcAlt,_pc,
        _knl,_pnl,_kpp,_pp,_kv0,_pv,_kvi,_pvi,_pixPar,_radStrength;
};

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Cosmology distorted power correlation function");
    std::string input,delta,output;
    int ellMax,nr,repeat,nk,nmu,samplesPerDecade,nrprt;
    double rmin,rmax,relerr,abserr,abspow,maxRelError,kmin,kmax,margin,vepsMin,vepsMax,drprt;
    double bias,biasbeta,biasGamma,biasSourceAbsorber,biasAbsorberResponse,meanFreePath,
        snlPar,snlPerp,kc,kcAlt,pc,knl,pnl,kpp,pp,kv0,pv,kvi,pvi,pixPar;
    cli.add_options()
        ("help,h", "prints this info and exits.")
        ("verbose", "prints additional information.")
        ("input,i", po::value<std::string>(&input)->default_value(""),
            "filename to read k,P(k) values from")
        ("delta", po::value<std::string>(&delta)->default_value(""),
            "optional filename of k,P(k) values to subtract from input")
        ("output,o", po::value<std::string>(&output)->default_value(""),
            "base name for saving results")
        ("rmin", po::value<double>(&rmin)->default_value(10.),
            "minimum value of comoving separation to use")
        ("rmax", po::value<double>(&rmax)->default_value(200.),
            "maximum value of comoving separation to use")
        ("nr", po::value<int>(&nr)->default_value(191),
            "number of points spanning [rmin,rmax] to use")
        ("ell-max", po::value<int>(&ellMax)->default_value(4),
            "maximum multipole to use for transforms")
        ("asymmetric", "distortion is asymmetric in mu (uses odd ell values)")
        ("bias", po::value<double>(&bias)->default_value(-0.17),
            "linear tracer bias")
        ("biasbeta", po::value<double>(&biasbeta)->default_value(-0.17),
            "product of bias and linear redshift-space distortion parameter beta")
        ("bias-gamma", po::value<double>(&biasGamma)->default_value(0),
            "Lyman-alpha photoionization bias")
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
        ("pix-par", po::value<double>(&pixPar)->default_value(0.),
            "scale for parallel pixelization in Mpc/h (ignored when pix-par = 0)")
        ("relerr", po::value<double>(&relerr)->default_value(1e-3),
            "relative error termination goal")
        ("abserr", po::value<double>(&abserr)->default_value(1e-5),
            "absolute error termination goal")
        ("abspow", po::value<double>(&abspow)->default_value(0.),
            "absolute error weighting power")
        ("margin", po::value<double>(&margin)->default_value(2.),
            "termination criteria margin to use for initialization")
        ("veps-max", po::value<double>(&vepsMax)->default_value(0.01),
            "maximum value of veps value to use")
        ("veps-min", po::value<double>(&vepsMin)->default_value(1e-6),
            "minimum value of veps value to use")
        ("samples-per-decade", po::value<int>(&samplesPerDecade)->default_value(40),
            "number of samples per decade to use for transform interpolation in k")
        ("max-rel-error", po::value<double>(&maxRelError)->default_value(1e-3),
            "maximum allowed relative error for power-law extrapolation of input P(k)")
        ("direct-power-multipoles",
            "use direct calculation of P(k) multipoles instead of interpolation")
        ("optimize", "optimizes transform FFTs")
        ("bypass", "bypasses the termination test for transforms")
        ("repeat", po::value<int>(&repeat)->default_value(1),
            "number of times to repeat identical transform")
        ("kmin", po::value<double>(&kmin)->default_value(0.005),
            "minimum value of comoving separation to use")
        ("kmax", po::value<double>(&kmax)->default_value(0.35),
            "maximum value of comoving separation to use")
        ("nk", po::value<int>(&nk)->default_value(100),
            "number of log-spaced k values for saving results")
        ("nmu", po::value<int>(&nmu)->default_value(10),
            "number of equally spaced mu_k and mu_r values for saving results")
        ("theta-angle", "use equally spaced theta angles for calculating mu_k and mu_r values")
        ("nrprt", po::value<int>(&nrprt)->default_value(0),
            "number of points along rp and rt axes for saving results")
        ("drprt", po::value<double>(&drprt)->default_value(4.),
            "spacing for points along rp and rt axes for saving results")
        ;
    // do the command line parsing now
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
    bool verbose(vm.count("verbose")), symmetric(0==vm.count("asymmetric")),
        optimize(vm.count("optimize")), bypass(vm.count("bypass")),
        directPowerMultipoles(vm.count("direct-power-multipoles")),
        thetaAngle(vm.count("theta-angle"));

    if(!symmetric) {
        std::cerr << "Odd multipoles not implemented yet." << std::endl;
        return 1;
    }

    if(input.length() == 0) {
        std::cerr << "Missing input filename." << std::endl;
        return 1;
    }

    int dell(symmetric ? 2:1);

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
            snlPar,snlPerp,kc,kcAlt,pc,knl,pnl,kpp,pp,kv0,pv,kvi,pvi,pixPar));
        cosmo::RMuFunctionCPtr distPtr(new cosmo::RMuFunction(boost::bind(
            &LyaDistortion::operator(),rsd,_1,_2)));

        // Use the limits of the input tabulated power for tabulating the
        // power multipoles (the kmin,kmax cmd-line args are for output only)
        double klo = power->getKMin(), khi = power->getKMax();
        int nkint = std::ceil(std::log10(khi/klo)*samplesPerDecade);
    	cosmo::DistortedPowerCorrelation dpc(PkPtr,distPtr,
            klo,khi,nkint,rmin,rmax,nr,ellMax,
            symmetric,relerr,abserr,abspow);
        // initialize
        dpc.initialize(nmu,margin,vepsMax,vepsMin,optimize);
        if(verbose) dpc.printToStream(std::cout);
        // transform (with repeats, if requested)
        bool ok;
        for(int i = 0; i < repeat; ++i) {
            ok = dpc.transform(!directPowerMultipoles,bypass);
        }
        if(!ok) {
            std::cerr << "Transform fails termination test." << std::endl;
        }
        // save transform results
        if(output.length() > 0) {
            int dell = symmetric ? 2 : 1;
            double mu;
            double dmu = symmetric ? 1./(nmu-1.) : 2./(nmu-1.);
            double dtheta = symmetric ? 2*std::atan(1)/(nmu-1.) : 4*std::atan(1)/(nmu-1.);
            // Write out values tabulated for log-spaced k
            double dk = std::pow(kmax/kmin,1./(nk-1.));
            std::string kfile = output + ".k.dat";
            std::ofstream kout(kfile.c_str());
            for(int i = 0; i < nk; ++i) {
                double k = kmin*std::pow(dk,i);
                kout << boost::lexical_cast<std::string>(k);
                for(int j = 0; j < nmu; ++j) {
                    if(thetaAngle) mu = std::cos(j*dtheta);
                    else mu = 1 - j*dmu;
                    kout << ' ' << boost::lexical_cast<std::string>(dpc.getPower(k,mu));
                }
                for(int ell = 0; ell <= ellMax; ell += dell) {
                    double pk = (directPowerMultipoles) ?
                        dpc.getPowerMultipole(k,ell) : dpc.getSavedPowerMultipole(k,ell);
                    kout << ' ' << boost::lexical_cast<std::string>(pk);
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
                    if(thetaAngle) mu = std::cos(j*dtheta);
                    else mu = 1 - j*dmu;
                    rout << ' ' << boost::lexical_cast<std::string>(dpc.getCorrelation(r,mu));
                }
                for(int ell = 0; ell <= ellMax; ell += dell) {
                    rout << ' ' << boost::lexical_cast<std::string>(dpc.getCorrelationMultipole(r,ell));
                }
                rout << std::endl;                
            }
            rout.close();
            // Write out values tabulated for (rp,rt) grid
            if(nrprt>0) {
                std::string rprtfile = output + ".rprt.dat";
                std::ofstream rprtout(rprtfile.c_str());
                for(int i = 0; i < nrprt; ++i) {
                    double rp = (0.5 + i)*drprt;
                    for(int j = 0; j < nrprt; ++j) {
                        double rt = (0.5 + j)*drprt;
                        double r = std::sqrt(rp*rp + rt*rt);
                        mu = rp/r;
                        rprtout << boost::lexical_cast<std::string>(rp) << ' ' <<
                        boost::lexical_cast<std::string>(rt) << ' ' << 
                        boost::lexical_cast<std::string>(dpc.getCorrelation(r,mu)) << std::endl;
                    }
                }
                rprtout.close();
            }
        }
    }
    catch(std::runtime_error const &e) {
        std::cerr << "ERROR: exiting with an exception:\n  " << e.what() << std::endl;
        return -1;
    }

    return 0;
}