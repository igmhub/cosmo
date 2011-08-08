# global compile and link options
AM_CPPFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = -llikely

# targets to build and install
lib_LTLIBRARIES = libcosmo.la
bin_PROGRAMS = cosmocalc

# extra targets that should not be installed
##noinst_PROGRAMS = demo1 demo2

# instructions for building the library
libcosmo_la_SOURCES = \
	cosmo/AbsHomogeneousUniverse.cc \
	cosmo/HomogeneousUniverseCalculator.cc \
	cosmo/LambdaCdmUniverse.cc

# library headers to install (nobase prefix preserves any subdirectories)
# Anything that includes config.h should *not* be listed here.
# and is not part of the package public API.
nobase_include_HEADERS = \
	cosmo/cosmo.h \
	cosmo/types.h \
	cosmo/RuntimeError.h \
	cosmo/AbsHomogeneousUniverse.h \
	cosmo/HomogeneousUniverseCalculator.h \
	cosmo/LambdaCdmUniverse.h

# instructions for building each program

cosmocalc_SOURCES = src/cosmocalc.cc
cosmocalc_DEPENDENCIES = $(lib_LIBRARIES)
cosmocalc_LDADD = -lboost_program_options -L. -lcosmo
