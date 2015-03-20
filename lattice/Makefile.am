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

check_PROGRAMS = graph_ck

graph_ck_SOURCES = graph_ck.cc

doc_DATA = lattice.pdf

SUFFIXES = .tex .pdf

EXTRA_DIST = lattice.tex

lattice.pdf: lattice.tex graph.tex

.hh.tex:
	../blib/glweb -o $@ -c $^

.tex.pdf:
	pdflatex $<

CLEANFILES = graph.tex *.aux *.log *.pdf

