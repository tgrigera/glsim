## 
## Makefile.am for test
##
## This file is part of glsim, a numerical simulation class library and
## helper programs.
##
## glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014
## by Tomas S. Grigera.

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

EXTRA_DIST += logging.nw
logging_hdr = 
logging_src = logtest.cc
BUILT_SOURCES += $(logging_hdr) $(logging_src)
$(logging_src) : logging.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(logging_hdr) : logging.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@

EXTRA_DIST += parameters.nw
parameters_hdr = 
parameters_src = partest.cc partest.ini partest2.ini
BUILT_SOURCES += $(parameters_hdr) $(parameters_src)
$(parameters_src) : parameters.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(parameters_hdr) : parameters.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@

EXTRA_DIST += environment.nw
environment_hdr = 
environment_src = envtest.cc envtest.ini envscope_test.cc envscope_test.ini
BUILT_SOURCES += $(environment_hdr) $(environment_src)
$(environment_src) : environment.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(environment_hdr) : environment.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@

EXTRA_DIST += random.nw
random_hdr = 
random_src = rantest.cc
BUILT_SOURCES += $(random_hdr) $(random_src)
$(random_src) : random.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(random_hdr) : random.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@

EXTRA_DIST += simulation.nw
simulation_hdr = 
simulation_src = simtest.cc simtest.ini
BUILT_SOURCES += $(simulation_hdr) $(simulation_src)
$(simulation_src) : simulation.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(simulation_hdr) : simulation.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@


EXTRA_DIST += $(BUILT_SOURCES)

AM_CPPFLAGS = -I$(top_builddir)/include
AM_CXXFLAGS=--std=c++11

check_PROGRAMS = context cerrors-test logtest partest envtest	\
	         envscope_test rantest simtest

context_SOURCES = context.cc ../blib/scontext.cc
cerrors_test_SOURCES = cerrors-test.c ../blib/cerrors.c
logtest_SOURCES = logtest.cc ../blib/log.cc
partest_SOURCES = partest.cc ../blib/parameters.cc ../blib/log.cc	\
	          ../blib/scontext.cc
envtest_SOURCES = envtest.cc ../blib/environment.cc	\
                  ../blib/parameters.cc ../blib/log.cc	\
                  ../blib/scontext.cc
envscope_test_SOURCES = envscope_test.cc ../blib/environment.cc	\
                  ../blib/parameters.cc ../blib/log.cc	\
                  ../blib/scontext.cc
rantest_SOURCES = rantest.cc ../blib/random.cc ../blib/scontext.cc
simtest_SOURCES = simtest.cc ../blib/simulation.cc		\
		  ../blib/environment.cc ../blib/parameters.cc	\
		  ../blib/log.cc ../blib/scontext.cc
