#
# Makefile.am - top-level automake input for glsim
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
# glsim is copyright (C) 2009, 2010, 2011, 2012, 2013 by Tomas S. Grigera.
# 
#


ACLOCAL_AMFLAGS=-I m4

.PHONY:  doc
if HASLATEX
doc:
	cd doc && make
else
doc:
	echo "***** LaTeX not found *****"
	echo "Cannot build documentation"
endif

idir = $(top_builddir)/include/glsim

$(idir):
	mkdir -p $@

BUILT_SOURCES = $(idir)

SUBDIRS = . blib test $(CONDITIONAL_SUBDIRS)

## SUBDIRS = base test $(CONDITIONAL_SUBDIRS) .
