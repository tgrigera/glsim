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

pkginclude_HEADERS = avevar.hh binvec.hh geoave.hh Grk.hh gr.hh	\
                     histogram.hh iscatt.hh Sk.hh timecorr.hh

EXTRA_DIST = gs_columnave.awk gs_fileave.awk

#
# library
#

lib_LTLIBRARIES = libglsim_an.la
libglsim_an_la_LDFLAGS = -version-info @LIBGLSIM_VERSION@
libglsim_an_la_SOURCES = avevar.cc geoave.cc gr.cc Grk.cc Sk.cc	\
                         timecorr.cc iscatt.cc

#
# binaries
#

bin_SCRIPTS = gs_columnave.awk gs_fileave.awk

bin_PROGRAMS = gs_1stneighbours gs_avesd gs_gav gs_gr gs_histogram	\
               gs_tcorr gs_tau gs_Fk gs_Sk

LDADD=libglsim_an.la $(top_builddir)/blib/libglsim.la

gs_1stneighbours_SOURCES = gs_1stneighbours.cc
gs_1stneighbours_LDADD = $(LDADD) $(top_builddir)/offlattice/libglsim_ol.la
gs_1stneighbours_CXXFLAGS=$(OPENMP_CXXFLAGS) $(AM_CXXFLAGS)

gs_avesd_SOURCES = avesd.cc

gs_Fk_SOURCES = Fk.cc
gs_Fk_LDADD=$(top_builddir)/offlattice/libglsim_ol.la libglsim_an.la	\
            $(top_builddir)/blib/libglsim.la

gs_gav_SOURCES = gav.cc

gs_gr_SOURCES = gs_gr.cc
gs_gr_LDADD=libglsim_an.la $(top_builddir)/offlattice/libglsim_ol.la \
            $(top_builddir)/blib/libglsim.la

gs_histogram_SOURCES = histocl.cc

gs_Sk_SOURCES = gs_Sk.cc
gs_Sk_LDADD=libglsim_an.la $(top_builddir)/offlattice/libglsim_ol.la \
            $(top_builddir)/blib/libglsim.la

gs_tcorr_SOURCES = tcorr.cc

gs_tau_SOURCES = tau.cc

