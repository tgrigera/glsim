## 
## Makefile.am for test
##
## This file is part of glsim, a numerical simulation class library and
## helper programs.
##
## glsim is copyright (C) 2009, 2010, 2011, 2012, 2013 by Tomas S. Grigera.

##idir = $(top_builddir)/offlattice/offlattice

EXTRA_DIST = test.nw

EXTRA_DIST += errors.nw
errors_hdr = 
errors_src = context.cc cerrors-test.c
BUILT_SOURCES = $(errors_hdr) $(errors_src)
$(errors_src) : errors.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(errors_hdr) : errors.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
	ln -s -f $(PWD)/$@ $(top_builddir)/include/glsim

EXTRA_DIST += logging.nw
logging_hdr = 
logging_src = logtest.cc
BUILT_SOURCES += $(logging_hdr) $(logging_src)
$(logging_src) : logging.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(logging_hdr) : logging.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
	ln -s -f $(PWD)/$@ $(top_builddir)/include/glsim

EXTRA_DIST += $(BUILT_SOURCES)

AM_CPPFLAGS = -I$(top_builddir)/include

check_PROGRAMS = context cerrors-test logtest

context_SOURCES = context.cc ../blib/scontext.cc
cerrors_test_SOURCES = cerrors-test.c ../blib/cerrors.c
logtest_SOURCES = logtest.cc ../blib/log.cc
