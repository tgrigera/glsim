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

## AM_CPPFLAGS = -I$(top_builddir)/include
AM_CXXFLAGS=--std=c++11

#
# Headers
#

pkginclude_HEADERS = cerrors.h exception.hh log.hh scontext.hh	\
                     random.hh configuration.hh environment.hh	\
                     simulation.hh observable.hh stochastic.hh

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
                      stochastic.cc observable.cc
