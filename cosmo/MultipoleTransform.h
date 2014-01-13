// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_MULTIPOLE_TRANSFORM
#define COSMO_MULTIPOLE_TRANSFORM

#include "boost/smart_ptr.hpp"

#include <vector>

namespace cosmo {
	class MultipoleTransform {
	// Calculates 2D (Hankel) or 3D (spherical Bessel) multipole transforms of
	// an arbitrary real-valued function.
	public:
		enum Type { SphericalBessel, Hankel };
		enum Strategy { EstimatePlan, MeasurePlan };
		// Creates a new transform object for an arbitrary func(u) that evaluates:
		//
		//   T(v) = Integrate[ S(ell,u,v)*func(u) , {u,0,Infinity} ]
		//
		// where S = u^2 j_ell(u*v) when type is SphericalBessel or
		// S = u J_ell(u*v) when type is Hankel. The result T(v) will be
		// evaluated to an accuracy over the range vmin < v < vmax that is
		// determined by the value of veps, with smaller values giving more
		// accurate results and requiring correspondingly more memory and cpu.
		// If veps > 0, then it roughly corresponds to 1/Nf, the number of
		// logarithmically spaced points where S is tabulated. Alternatively,
		// when veps < 0 then the symmetrized S' is truncated at smax with
		// S'(smax) = (-veps)*S'(0).
		MultipoleTransform(Type type, int ell, double vmin, double vmax, double veps,
			Strategy strategy, int minSamplesPerCycle = 2, int minSamplesPerDecade = 40);
		// Creates a new transform object that has parameters as an existing transform
		// object except that it uses a logarithmic spacing that is larger by the
		// specified integer subsampling factor.
		MultipoleTransform(MultipoleTransform const &other, int subsampling);
		virtual ~MultipoleTransform();
		// Returns the truncation fraction eps such that the symmetrized S' is
		// assumed to be zero for |s| > smax with S'(smax) = eps*S'(0). This is the
		// same value that can specified directly in the constructor using veps = -eps.
		double getTruncationFraction() const;
		// Returns the minimum number of samples per cycle for this transformer.
		int getMinSamplesPerCycle() const;
		// Returns the grid of u values where a function to be transformed will be
		// evaluated when calling the transform(...) method.
		std::vector<double> const &getUGrid() const;
		// Returns the grid of v values where our transformed result will be estimated
		// after calling the transform(...) method. Note that the range of v values will
		// generally be larger than [vmin,vmax].
		std::vector<double> const &getVGrid() const;
		// Estimates the transform of func on our v grid using the the specified
		// values of func(u) tabulated on our u grid. The results are saved in
		// the results vector provided, which will be resized if necesary.
		void transform(std::vector<double> const &funcTable,
			std::vector<double> &result) const;
	private:
		Type _type;
		double _eps;
		int _minSamplesPerCycle;
		std::vector<double> _ugrid, _vgrid, _coef, _scale;
		// We use an implementation subclass to avoid any public include dependency
		// on fftw, since this is an optional package when building our library.
		class Implementation;
		boost::scoped_ptr<Implementation> _pimpl;
	}; // MultipoleTransform

	inline double MultipoleTransform::getTruncationFraction() const {
		return _eps;
	}
	inline int MultipoleTransform::getMinSamplesPerCycle() const {
		return _minSamplesPerCycle;
	}
	inline std::vector<double> const &MultipoleTransform::getUGrid() const {
		return _ugrid;
	}
	inline std::vector<double> const &MultipoleTransform::getVGrid() const {
		return _vgrid;
	}

} // cosmo

#endif // COSMO_MULTIPOLE_TRANSFORM
