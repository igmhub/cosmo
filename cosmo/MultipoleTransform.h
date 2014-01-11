// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_MULTIPOLE_TRANSFORM
#define COSMO_MULTIPOLE_TRANSFORM

#include "likely/function.h"

#include <vector>

namespace cosmo {
	class MultipoleTransform {
	// Calculates 2D (Hankel) and 3D (spherical Bessel) multipole transforms of
	// an arbitrary real-valued function.
	public:
		enum Type { SphericalBessel, Hankel };
		// Creates a new transform object for an arbitrary func(u) that evaluates:
		//
		//   T(v) = Integrate[ S(ell,u,v)*func(u) , {u,0,Infinity} ]
		//
		// where S = u^2 j_ell(u*v) when type is SphericalBessel or
		// S = u J_ell(u*v) when type is Hankel. The result T(v) will be
		// evaluated to the accuracy eps over the range vmin < v < vmax.
		MultipoleTransform(Type type, int ell,
			double vmin, double vmax, double veps, int minSamplesPerDecade = 40);
		virtual ~MultipoleTransform();
		// Calculates the transform of the specified function on the specified
		// v grid, saving the results in the vector provided (which will be
		// resized if necessary).
		void transform(likely::GenericFunctionPtr func, std::vector<double> const &vgrid,
			std::vector<double> &result) const;
	private:
	}; // MultipoleTransform
} // cosmo

#endif // COSMO_MULTIPOLE_TRANSFORM
