#
# Makefile.am for lattice directory
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
AM_CPPFLAGS = -I$(top_srcdir)/lib/blib
LDADD=$(top_builddir)/lib/libglsim.la

#
# Headers
#

pkginclude_HEADERS = graph.hh fcgraph.hh bethe.hh lattice1D.hh	\
                     lattice2D.hh lattice3D.hh graph_test.hh
if HAS_HDF
pkginclude_HEADERS += poisson.hh
endif

#
# binaries
#
check_PROGRAMS = graph_ck fcgraph_test bethe_test lattice1D_test lattice2D_test	\
                 lattice3D_test
if HAS_HDF
check_PROGRAMS += poisson_test
endif

# Binary dependencies

graph_ck_SOURCES = graph_ck.cc
fcgraph_test_SOURCES = graph_test.cc fcgraph_test.cc
bethe_test_SOURCES = graph_test.cc bethe_test.cc
lattice1D_test_SOURCES = graph_test.cc lattice1D_test.cc
lattice2D_test_SOURCES = graph_test.cc lattice2D_test.cc
lattice3D_test_SOURCES = graph_test.cc lattice3D_test.cc

if HAS_HDF
poisson_test_SOURCES = graph_test.cc poisson_test.cc
poisson_test_CPPFLAGS = -I$(top_srcdir)/lib/blib -I$(top_srcdir)/lib/offlattice
poisson_test_LDADD = $(top_builddir)/lib/libglsim.la $(top_builddir)/lib/libglsim_ol.la
endif


#
# Tests
#

TESTS = fcgraph_test bethe_test lattice1D_test lattice2D_test lattice3D_test

if HAS_HDF
TESTS += poisson_test
endif


