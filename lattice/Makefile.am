#
# Makefile.am for lattice directory
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
LDADD=$(top_builddir)/blib/libglsim.la

#
# Headers
#

pkginclude_HEADERS = graph.hh bethe.hh lattice1D.hh lattice2D.hh	\
                     lattice3D.hh test_exception.hh graph_test.hh

#
# binaries
#
check_PROGRAMS = graph_ck bethe_test lattice1D_test lattice2D_test	\
                 lattice3D_test

# Binary dependencies

graph_ck_SOURCES = graph_ck.cc
bethe_test_SOURCES = graph_test.cc bethe_test.cc
lattice1D_test_SOURCES = graph_test.cc lattice1D_test.cc
lattice2D_test_SOURCES = graph_test.cc lattice2D_test.cc
lattice3D_test_SOURCES = graph_test.cc lattice3D_test.cc

#
# Tests
#

TESTS = bethe_test lattice1D_test lattice2D_test lattice3D_test

#doc_DATA = lattice.pdf

#CLEANFILES = graph.tex *.aux *.log *.pdf

