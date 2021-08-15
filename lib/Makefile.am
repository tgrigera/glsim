#
# Makefile.am for blib directory
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
#  glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
#                         2017, 2018, 2019
#  by Tomas S. Grigera.
# 
# glsim is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License (GPL) as published
# by the Free Software Foundation. You can use either version 3, or
# (at your option) any later version.
#
# If you use glsim to produced published work, or if you redistribute
# a modified version of glsim, or code based on glsim, please cite us
#
# - T. S. Grigera, /glsim: A general library for numerical
#   simulation/.  Computer Physics Communications *182*, 2122-2131
#   (2011).
#
# glsim is developed as part of academic work.  As such author
# attributions are not only a way to give recognition to the original
# authors, but also help to ensure continued development, since
# citations show that the work has been useful to others and thus
# help convince funding agencies that it is worth funding.  * glsim
# distribution.
#
# glsim is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.
# 
# For details see the file LICENSE in the home directory. If the file is
# missing, contact the maintainers.
#

#
# General compiler and linker options
#

include $(top_srcdir)/Makefile-common
AUTOMAKE_OPTIONS = subdir-objects
AM_CPPFLAGS = -I$(top_srcdir)/lib/blib

#
# Basic library
#

# Headers

pkginclude_HEADERS = blib/blib.hh blib/cerrors.h blib/exception.hh	\
                     blib/fft.hh blib/log.hh blib/parameters.hh		\
                     blib/scontext.hh blib/random.hh			\
                     blib/configuration.hh blib/environment.hh		\
                     blib/simulation.hh blib/observable.hh		\
                     blib/stochastic.hh blib/mfile.hh			\
                     blib/test_exception.hh

pkginclude_HEADERS += analysis/avevar.hh analysis/binvec.hh	\
                     analysis/geoave.hh 		\
                     analysis/histogram.hh analysis/iscatt.hh	\
                     analysis/timecorr.hh


# Doc-only headers

EXTRA_DIST = blib-doc.hh

# library

lib_LTLIBRARIES = libglsim.la
libglsim_la_LDFLAGS = -version-info @LIBGLSIM_VERSION@
libglsim_la_SOURCES = blib/scontext.cc blib/cerrors.c blib/log.cc	\
                      blib/random.cc blib/parameters.cc			\
                      blib/environment.cc blib/simulation.cc		\
                      blib/stochastic.cc blib/observable.cc		\
                      blib/mfile.cc blib/fft.cc

# analysis modules
libglsim_la_SOURCES += analysis/avevar.cc analysis/geoave.cc	\
                       analysis/timecorr.cc analysis/iscatt.cc

# Test programs

LDADD=libglsim.la

check_PROGRAMS = exception_test parameters_test environment_test	\
                 fft_test mvector_test

exception_test_SOURCES = exception_test.cc

environment_test_SOURCES = environment_test.cc

fft_test_SOURCES = fft_test.cc

parameters_test_SOURCES = parameters_test.cc

mvector_test_SOURCES = mvector_test.cc

# Tests

TESTS = fft_test


#
# Offlattice lib
#

if HAS_HDF

AM_CPPFLAGS += -I$(top_srcdir)/lib/offlattice

LDADD += libglsim_ol.la

# offlattice headers

pkginclude_HEADERS += offlattice/hdf_file.hh offlattice/offlattice.hh	\
		     offlattice/h5md.hh offlattice/olconfiguration.hh	\
		     offlattice/interactions.hh offlattice/ld.hh	\
		     offlattice/lj.hh offlattice/mc.hh			\
		     offlattice/mcobservable.hh				\
		     offlattice/mdenvironment.hh			\
		     offlattice/mdobservable.hh offlattice/md.hh	\
		     offlattice/nneighbours.hh				\
		     offlattice/trajectory.hh

pkginclude_HEADERS += analysis/Grk.hh analysis/Sk.hh analysis/gr.hh

# offlattice lib 

lib_LTLIBRARIES += libglsim_ol.la
libglsim_ol_la_LDFLAGS = -version-info @LIBGLSIM_VERSION@
libglsim_ol_la_SOURCES = offlattice/hdf_file.cc offlattice/h5md.cc	\
			 offlattice/olconfiguration.cc			\
			 offlattice/interactions.cc offlattice/ld.cc	\
			 offlattice/mc.cc offlattice/mcobservable.cc	\
			 offlattice/md.cc offlattice/mdenvironment.cc	\
			 offlattice/mdobservable.cc			\
			 offlattice/nneighbours.cc			\
			 offlattice/trajectory.cc

# offlattice analysis

libglsim_ol_la_SOURCES += analysis/Grk.cc analysis/Sk.cc analysis/gr.cc

# offlattice tests

check_PROGRAMS += hdf_recfile_test neighbour_test pseudo-md

hdf_recfile_test_SOURCES = hdf_recfile_test.cc
neighbour_test_SOURCES = neighbour_test.cc
pseudo_md_SOURCES = pseudo-md.cc

TESTS += neighbour_test


endif
