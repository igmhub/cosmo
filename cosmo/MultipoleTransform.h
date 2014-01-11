// Created 10-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_MULTIPOLE_TRANSFORM
#define COSMO_MULTIPOLE_TRANSFORM

#include "likely/function.h"

namespace cosmo {
	class MultipoleTransform {
	// Calculates 2D (Hankel) and 3D (spherical Bessel) multipole transforms of
	// an arbitrary real-valued function.
	public:
		enum Type { SphericalBessel, Hankel };
		// Creates a new transform object for the specified func(u) that evaluates:
		//
		//   T(v) = Integrate[ S(ell,u,v)*func(u) , {u,0,Infinity} ]
		//
		// where S = u^2 j_ell(u*v) when type is SphericalBessel or
		// S = u J_ell(u*v) when type is Hankel. The result T(v) will be
		// evaluated to the accuracy eps over the range vmin < v < vmax.
		MultipoleTransform(likely::GenericFunctionPtr func, Type type, int ell,
			double vmin, double vmax, double veps, int minSamplesPerDecade = 40);
		virtual ~MultipoleTransform();
	private:
	}; // MultipoleTransform
} // cosmo

#endif // COSMO_MULTIPOLE_TRANSFORM
