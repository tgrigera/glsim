## 
## Makefile.am for test
##
## This file is part of glsim, a numerical simulation class library and
## helper programs.
##
## glsim is copyright (C) 2009, 2010, 2011, 2012, 2013 by Tomas S. Grigera.

errors_hdr = scontext.hh exception.hh cerrors.h
errors_src = scontext.cc cerrors.c
errors_scrpt = readbt.sh
BUILT_SOURCES = $(errors_hdr) $(errors_src) $(errors_scrpt)
EXTRA_DIST = errors.nw $(BUILT_SOURCES)
$(errors_src) : errors.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(errors_hdr) : errors.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
	ln -s -f $(PWD)/$@ $(top_builddir)/include/glsim
$(errors_scrpt) : errors.nw
	@TANGLE@ -R$@ $? | @CPIF@ $@
	chmod +x $@

pkginclude_HEADERS = $(errors_hdr)

bin_SCRIPTS = readbt.sh
