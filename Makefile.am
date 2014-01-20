# add path to boost.m4
ACLOCAL_AMFLAGS = -I .

# global compile and link options
AM_CPPFLAGS = $(BOOST_CPPFLAGS)

# targets to build and install
lib_LTLIBRARIES = libcosmo.la
bin_PROGRAMS = cosmocalc cosmo3d cosmogrf cosmostack cosmoxi cosmomock \
cosmotrans cosmoatrans

# extra targets that should not be installed
noinst_PROGRAMS = cosmotest

# targets that contain unit tests
#check_PROGRAMS = 

# add our pkgconfig file to the install target
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = cosmo.pc

# any library dependencies not already added by configure can be added here
#libcosmo_la_LIBADD = 

# instructions for building the library
libcosmo_la_SOURCES = \
	cosmo/AbsHomogeneousUniverse.cc \
	cosmo/HomogeneousUniverseCalculator.cc \
	cosmo/LambdaCdmUniverse.cc \
	cosmo/LambdaCdmRadiationUniverse.cc \
	cosmo/BaryonPerturbations.cc \
	cosmo/BroadbandPower.cc \
	cosmo/TabulatedPower.cc \
	cosmo/TransferFunctionPowerSpectrum.cc \
	cosmo/PowerSpectrumCorrelationFunction.cc \
	cosmo/RsdCorrelationFunction.cc \
	cosmo/OneDimensionalPowerSpectrum.cc \
	cosmo/AbsGaussianRandomFieldGenerator.cc \
	cosmo/FftGaussianRandomFieldGenerator.cc \
	cosmo/TestFftGaussianRandomFieldGenerator.cc \
	cosmo/MultipoleTransform.cc \
	cosmo/AdaptiveMultipoleTransform.cc \
	cosmo/DistortedPowerCorrelation.cc

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
	cosmo/TabulatedPower.h \
	cosmo/TransferFunctionPowerSpectrum.h \
	cosmo/PowerSpectrumCorrelationFunction.h \
	cosmo/RsdCorrelationFunction.h \
	cosmo/OneDimensionalPowerSpectrum.h \
	cosmo/AbsGaussianRandomFieldGenerator.h \
	cosmo/FftGaussianRandomFieldGenerator.h \
	cosmo/TestFftGaussianRandomFieldGenerator.h \
	cosmo/MultipoleTransform.h \
	cosmo/AdaptiveMultipoleTransform.h \
	cosmo/DistortedPowerCorrelation.h

# instructions for building each program

cosmocalc_SOURCES = src/cosmocalc.cc
cosmocalc_DEPENDENCIES = $(lib_LIBRARIES)
cosmocalc_LDADD = libcosmo.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)

cosmo3d_SOURCES = src/cosmo3d.cc
cosmo3d_DEPENDENCIES = $(lib_LIBRARIES)
cosmo3d_LDADD = libcosmo.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)

cosmotest_SOURCES = src/cosmotest.cc
cosmotest_DEPENDENCIES = $(lib_LIBRARIES)
cosmotest_LDADD = libcosmo.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)

cosmogrf_SOURCES = src/cosmogrf.cc
cosmogrf_DEPENDENCIES = $(lib_LIBRARIES)
cosmogrf_LDADD = libcosmo.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)

cosmostack_SOURCES = src/cosmostack.cc
cosmostack_DEPENDENCIES = $(lib_LIBRARIES)
cosmostack_LDADD = libcosmo.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)

cosmoxi_SOURCES = src/cosmoxi.cc
cosmoxi_DEPENDENCIES = $(lib_LIBRARIES)
cosmoxi_LDADD = libcosmo.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)

cosmomock_SOURCES = src/cosmomock.cc
cosmomock_DEPENDENCIES = $(lib_LIBRARIES)
cosmomock_LDADD = libcosmo.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)

cosmotrans_SOURCES = src/cosmotrans.cc
cosmotrans_DEPENDENCIES = $(lib_LIBRARIES)
cosmotrans_LDADD = libcosmo.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)

cosmoatrans_SOURCES = src/cosmoatrans.cc
cosmoatrans_DEPENDENCIES = $(lib_LIBRARIES)
cosmoatrans_LDADD = libcosmo.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
