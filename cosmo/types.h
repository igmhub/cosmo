// Created 8-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_TYPES
#define COSMO_TYPES

#include "boost/smart_ptr.hpp"

namespace cosmo {
    
    class AbsHomogeneousUniverse;
    typedef boost::shared_ptr<AbsHomogeneousUniverse> AbsHomogeneousUniversePtr;
    
} // cosmo

#endif COSMO_TYPES
