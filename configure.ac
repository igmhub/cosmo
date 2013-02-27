#################################################################################
## Created 8-Aug-2011 by David Kirkby <dkirkby@uci.edu>
## Run autoreconf after making any change here.
#################################################################################

# This could probably be relaxed to an older version.
AC_PREREQ([2.68])

# Our package info.
AC_INIT([cosmo], [dev], [dkirkby@uci.edu])

# Checks for programs
AC_PROG_CXX

# Initialize libtool, which adds --enable/disable-shared configure options.
# The boost.m4 macros used below also need this.
LT_INIT

# Checks for libraries (there is no portable way to check for C++ classes in
# the library, so we just check that the linker can find the library using 'main')
AC_CHECK_LIB([m],[cos],,
	AC_MSG_ERROR([Cannot find the math library.]))

# Use 'configure --without-fftw3' if no FFTW3 library is available.
AC_ARG_WITH([fftw3],
	AS_HELP_STRING([--without-fftw3], [Build without the FFTW3 library.]))
AS_IF([test "x$with_fftw3" != "xno"], [
	AC_CHECK_LIB([fftw3f],[fftwf_malloc],,
		AC_MSG_ERROR([Cannot find the FFTW3 single-precision library.]))
])

# We need a recent version of boost
BOOST_REQUIRE([1.49])

# Required header-only boost packages
BOOST_BIND
BOOST_FOREACH
BOOST_FUNCTION
BOOST_UTILITY
BOOST_SMART_PTR

# Required boost libraries
BOOST_PROGRAM_OPTIONS

# Configure automake
AC_CONFIG_FILES([Makefile])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
# See http://www.gnu.org/s/hello/manual/automake/maintainer_002dmode.html
AM_MAINTAINER_MODE
AC_CONFIG_HEADERS([config.h])

# Define automake variables that flag whether optional libraries should be used.
AM_CONDITIONAL([USE_FFTW3], [test "x$ac_cv_lib_fftw3_fftw_malloc" = "xyes"])

AC_OUTPUT
