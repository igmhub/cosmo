// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_MULTIPOLE_TRANSFORM
#define COSMO_MULTIPOLE_TRANSFORM

#include "boost/smart_ptr.hpp"

#include <vector>

namespace cosmo {
	class MultipoleTransform {
	// Calculates 2D (Hankel) or 3D (spherical Bessel) multipole transforms of
	// an arbitrary real-valued function. For details on the method, see
	// https://www.authorea.com/users/4112/articles/4271
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
		// determined by the value of |veps|, with smaller values giving more
		// accurate results and requiring correspondingly more memory and cpu.
		// If veps > 0, then it roughly corresponds to 1/Nf, the number of
		// logarithmically spaced points where S is tabulated. Alternatively,
		// when veps < 0 then the symmetrized S' is truncated at smax with
		// S'(smax) = (-veps)*S'(0). The strategy selects a tradeoff between
		// initialization and transform speeds (via the FFTW plan strategy option).
		// Different strategies can give different numerical results at the level
		// of roundoff errors.
		MultipoleTransform(Type type, int ell, double vmin, double vmax, double veps,
			Strategy strategy, int minSamplesPerCycle = 2, int minSamplesPerDecade = 40,
			int interpolationPadding = 3);
		virtual ~MultipoleTransform();
		// Returns the truncation fraction eps such that the symmetrized S' is
		// assumed to be zero for |s| > smax with S'(smax) = eps*S'(0). This is the
		// same value that can specified directly in the constructor using veps = -eps.
		double getTruncationFraction() const;
		// Returns the minimum number of samples per cycle for this transformer.
		int getMinSamplesPerCycle() const;
		// Returns the number of logarithmically spaced points where the symmetrized S'
		// is evaluated for convolution. Note that this is less than the size of our
		// u grid because of the zero padding that is added to eliminate aliasing artifacts.
		int getNumPoints() const;
		// Returns the number of logarithmically-spaced samples per decade used to
		// evaluate the symmetrized S' for convolution.
		double getSamplesPerDecade() const;
		// Returns the grid of u values where a function to be transformed should be
		// evaluated when preparing the funcTable for calling the transform(...) method.
		// Note that the values in u grid are always strictly decreasing (!) and positive.
		std::vector<double> const &getUGrid() const;
		// Returns the grid of v values where our transformed result will be estimated
		// after calling the transform(...) method. The algorithm internally uses a
		// range of v values that is much larger than [vmin,vmax] but this function only
		// returns the subrange that is guaranteed to be free of convolution aliasing
		// artifacts, and also guaranteed to extend beyond [vmin,vmax] by at least
		// interpolationPadding points on each side.
		std::vector<double> const &getVGrid() const;
		// Estimates the transform of func on our v grid using the the specified
		// values of func(u) tabulated on our u grid. The results are saved in
		// the results vector provided, which will be resized to our vgrid size
		// if necessary.
		void transform(std::vector<double> const &funcTable,
			std::vector<double> &result) const;
	private:
		Type _type;
		double _eps;
		int _minSamplesPerCycle, _Nf, _cleanBegin, _cleanEnd;
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
	inline int MultipoleTransform::getNumPoints() const {
		return 2*_Nf;
	}
	inline std::vector<double> const &MultipoleTransform::getUGrid() const {
		return _ugrid;
	}
	inline std::vector<double> const &MultipoleTransform::getVGrid() const {
		return _vgrid;
	}

	// Returns the coefficient of the Hankel (ndim = 2) or spherical Bessel (ndim = 3)
	// transform of...
	//   ...f_ell(r) that gives f~_ell(k) [dir = +1]
	//   ...f~_ell(k) that gives f_ell(r) [dir = -1]
	//
	// The assumed Fourier convention is specified by the parameters (a,b) such that:
	//
	//   f(vec{r}) = N_n(a,b)   Integral[ d^nk exp(+b(k.r)) f~(vec{k}) ]
	//   f~(vec{k}) = N~_n(a,b) Integral[ d^nr exp(+b(k.r)) f(vec{r}) ]
	//
	// where:
	//
	//   N_n(a,b)  = |b|^(n/2) (2pi)^(-n(1+a)/2)
	//   N~_n(a,b) = |b|^(/n2) (2pi)^(-n(1-a)/2)
	//
	// If ell is even, then the coefficient is real. However, if ell is odd, then the
	// coefficient is pure imaginary and the return value should be multiplied by i.
	double multipoleTransformNormalization(int ell, int ndim, int dir,
		double a = 1, double b = 1);

} // cosmo

#endif // COSMO_MULTIPOLE_TRANSFORM
