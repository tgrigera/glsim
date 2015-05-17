#
# Makefile.am - top-level automake input for glsim
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
# glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015
#                        by Tomas S. Grigera.
# 
#


ACLOCAL_AMFLAGS=-I m4

.PHONY:  doc
if HAS_DOXYGEN

CONDITIONAL_SUBDIRS = doc

doc:
	cd doc && make
else
doc:
	echo "***** Doxygen not found *****"
	echo "Cannot build documentation"
endif

# idir = $(top_builddir)/include/glsim

# $(idir):
# 	mkdir -p $@

# BUILT_SOURCES = $(idir)

SUBDIRS = . blib lattice $(CONDITIONAL_SUBDIRS)

MAINTAINERCLEANFILES = config.guess config.sub configure depcomp	\
		       install-sh ltmain.sh Makefile.in missing aclocal.m4 m4/*
