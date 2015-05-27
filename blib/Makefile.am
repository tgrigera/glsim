#
# Makefile.am for blib directory
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
# glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015
# by Tomas S. Grigera.
#
# See file LICENSE in base directory
#

#
# General compiler and linker options
#

include $(top_srcdir)/Makefile-common

#
# Headers
#

pkginclude_HEADERS = cerrors.h exception.hh log.hh parameters.hh	\
                     scontext.hh random.hh configuration.hh		\
                     environment.hh simulation.hh observable.hh		\
                     stochastic.hh

#
# Doc-only headers
#

EXTRA_DIST = blib.hh readbt.sh

#
# binaries
#

bin_SCRIPTS = readbt.sh

#
# library
#

lib_LTLIBRARIES = libglsim.la
libglsim_la_LDFLAGS = -version-info @LIBGLSIM_VERSION@
libglsim_la_SOURCES = scontext.cc cerrors.c log.cc random.cc		\
                      parameters.cc environment.cc simulation.cc	\
                      stochastic.cc observable.cc hdf_file.cc

#
# Test programs
#

LDADD=libglsim.la

check_PROGRAMS = exception_test parameters_test environment_test hdf_recfile_test

exception_test_SOURCES = exception_test.cc

environment_test_SOURCES = environment_test.cc

hdf_recfile_test_SOURCES = hdf_recfile_test.cc

parameters_test_SOURCES = parameters_test.cc
