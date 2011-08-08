// Created 08-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef COSMO_RUNTIME_ERROR
#define COSMO_RUNTIME_ERROR

#include <stdexcept>
#include <string>

namespace cosmo {
	class RuntimeError : public std::runtime_error {
	public:
		explicit RuntimeError(std::string const &reason);
		virtual ~RuntimeError() throw ();
	private:
	}; // RuntimeError
	
	inline RuntimeError::RuntimeError(std::string const &reason)
	: std::runtime_error(reason) { }

    inline RuntimeError::~RuntimeError() throw () { }
} // cosmo

#endif // COSMO_RUNTIME_ERROR
