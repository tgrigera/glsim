#
# Makefile.am for analysis directory
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
AM_CPPFLAGS = -I$(top_srcdir)/blib -I$(top_srcdir)/offlattice

#
# Headers
#

pkginclude_HEADERS = avevar.hh geoave.hh histogram.hh Sk.hh timecorr.hh

#
# library
#

lib_LTLIBRARIES = libglsim_an.la
libglsim_an_la_LDFLAGS = -version-info @LIBGLSIM_VERSION@
libglsim_an_la_SOURCES = avevar.cc geoave.cc Sk.cc timecorr.cc

#
# binaries
#

bin_SCRIPTS = gs_columnave.awk gs_fileave.awk

bin_PROGRAMS = gs_avesd gs_gav gs_histogram gs_tcorr

LDADD=libglsim_an.la $(top_builddir)/blib/libglsim.la

gs_avesd_SOURCES = avesd.cc

gs_gav_SOURCES = gav.cc

gs_histogram_SOURCES = histocl.cc

gs_tcorr_SOURCES = tcorr.cc
