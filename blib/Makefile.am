## 
## Makefile.am for test
##
## This file is part of glsim, a numerical simulation class library and
## helper programs.
##
## glsim is copyright (C) 2009, 2010, 2011, 2012, 2013 by Tomas S. Grigera.

EXTRA_DIST = blib.nw

EXTRA_DIST += errors.nw
errors_hdr = scontext.hh exception.hh cerrors.h
errors_src = scontext.cc cerrors.c
errors_scrpt = readbt.sh
BUILT_SOURCES = $(errors_hdr) $(errors_src) $(errors_scrpt)
$(errors_src) : errors.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(errors_hdr) : errors.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
	ln -s -f $(PWD)/$@ $(top_builddir)/include/glsim
$(errors_scrpt) : errors.nw
	@TANGLE@ -R$@ $? | @CPIF@ $@
	chmod +x $@

EXTRA_DIST += logging.nw
logging_hdr = log.hh
logging_src = log.cc
BUILT_SOURCES += $(logging_hdr) $(logging_src)
$(logging_src) : logging.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(logging_hdr) : logging.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
	ln -s -f $(PWD)/$@ $(top_builddir)/include/glsim

EXTRA_DIST += parameters.nw
parameters_hdr = parameters.hh
parameters_src = parameters.cc
BUILT_SOURCES += $(parameters_hdr) $(parameters_src)
$(parameters_src) : parameters.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(parameters_hdr) : parameters.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
	ln -s -f $(PWD)/$@ $(top_builddir)/include/glsim

EXTRA_DIST += $(BUILT_SOURCES)

pkginclude_HEADERS = $(errors_hdr) $(logging_hdr) $(parameters_hdr)

bin_SCRIPTS = readbt.sh
