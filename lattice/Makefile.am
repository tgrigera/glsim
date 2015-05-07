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

AM_CPPFLAGS = -I$(top_builddir)/include
AM_CXXFLAGS=--std=c++11

pkginclude_HEADERS = graph.hh

check_PROGRAMS = graph_ck lattice2D_test lattice3D_test

graph_ck_SOURCES = graph_ck.cc
lattice2D_test_SOURCES = graph_test.cc lattice2D_test.cc ../blib/scontext.cc
lattice3D_test_SOURCES = graph_test.cc lattice3D_test.cc ../blib/scontext.cc

TESTS = lattice2D_test lattice3D_test

#doc_DATA = lattice.pdf

#CLEANFILES = graph.tex *.aux *.log *.pdf

