# global compile and link options
AM_CPPFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = -llikely

# targets to build and install
lib_LTLIBRARIES = libcosmo.la
bin_PROGRAMS = cosmocalc cosmo3d

# extra targets that should not be installed
noinst_PROGRAMS = cosmotest

# instructions for building the library
libcosmo_la_SOURCES = \
	cosmo/AbsHomogeneousUniverse.cc \
	cosmo/HomogeneousUniverseCalculator.cc \
	cosmo/LambdaCdmUniverse.cc \
	cosmo/LambdaCdmRadiationUniverse.cc \
	cosmo/BaryonPerturbations.cc \
	cosmo/BroadbandPower.cc \
	cosmo/TransferFunctionPowerSpectrum.cc \
	cosmo/PowerSpectrumCorrelationFunction.cc \
	cosmo/RsdCorrelationFunction.cc \
	cosmo/OneDimensionalPowerSpectrum.cc \
	cosmo/AbsGaussianRandomFieldGenerator.cc \
	cosmo/FftGaussianRandomFieldGenerator.cc

# library headers to install (nobase prefix preserves any subdirectories)
# Anything that includes config.h should *not* be listed here.
# and is not part of the package public API.
nobase_include_HEADERS = \
	cosmo/cosmo.h \
	cosmo/types.h \
	cosmo/RuntimeError.h \
	cosmo/AbsHomogeneousUniverse.h \
	cosmo/HomogeneousUniverseCalculator.h \
	cosmo/LambdaCdmUniverse.h \
	cosmo/LambdaCdmRadiationUniverse.h \
	cosmo/BaryonPerturbations.h \
	cosmo/BroadbandPower.h \
	cosmo/TransferFunctionPowerSpectrum.h \
	cosmo/PowerSpectrumCorrelationFunction.h \
	cosmo/RsdCorrelationFunction.h \
	cosmo/OneDimensionalPowerSpectrum.h \
	cosmo/AbsGaussianRandomFieldGenerator.h \
	cosmo/FftGaussianRandomFieldGenerator.h

# instructions for building each program

cosmocalc_SOURCES = src/cosmocalc.cc
cosmocalc_DEPENDENCIES = $(lib_LIBRARIES)
cosmocalc_LDADD = -lboost_program_options -L. -lcosmo

cosmo3d_SOURCES = src/cosmo3d.cc
cosmo3d_DEPENDENCIES = $(lib_LIBRARIES)
cosmo3d_LDADD = -lboost_program_options -L. -lcosmo

cosmotest_SOURCES = src/cosmotest.cc
cosmotest_DEPENDENCIES = $(lib_LIBRARIES)
cosmotest_LDADD = -lboost_program_options -L. -lcosmo
