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
AM_CPPFLAGS = -I$(top_srcdir)/blib

#
# binaries
#

bin_SCRIPTS = gs_columnave.awk

bin_PROGRAMS = gs_avesd

LDADD=$(top_builddir)/blib/libglsim.la

gs_avesd_SOURCES = avesd.cc
