## 
## Makefile.am for test
##
## This file is part of glsim, a numerical simulation class library and
## helper programs.
##
## glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014
## by Tomas S. Grigera.

#
# Build sources from noweb files
#

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

EXTRA_DIST += random.nw
random_hdr = random.hh
random_src = random.cc
BUILT_SOURCES += $(random_hdr) $(random_src)
$(random_src) : random.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(random_hdr) : random.nw
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

EXTRA_DIST += environment.nw
environment_hdr = environment.hh ctenvironment.hh
environment_src = environment.cc ctenvironment.cc
BUILT_SOURCES += $(environment_hdr) $(environment_src)
$(environment_src) : environment.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(environment_hdr) : environment.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
	ln -s -f $(PWD)/$@ $(top_builddir)/include/glsim

EXTRA_DIST += configuration.nw
configuration_hdr = configuration.hh
configuration_src = configuration.cc
BUILT_SOURCES += $(configuration_hdr) $(configuration_src)
$(configuration_src) : configuration.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(configuration_hdr) : configuration.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
	ln -s -f $(PWD)/$@ $(top_builddir)/include/glsim

EXTRA_DIST += simulation.nw
simulation_hdr = simulation.hh
simulation_src = simulation.cc
BUILT_SOURCES += $(simulation_hdr) $(simulation_src)
$(simulation_src) : simulation.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(simulation_hdr) : simulation.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
	ln -s -f $(PWD)/$@ $(top_builddir)/include/glsim

EXTRA_DIST += stochastic.nw
stochastic_hdr = stochastic.hh
stochastic_src = stochastic.cc
BUILT_SOURCES += $(stochastic_hdr) $(stochastic_src)
$(stochastic_src) : stochastic.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
$(stochastic_hdr) : stochastic.nw
	@TANGLE@ -L -R$@ $? | @CPIF@ $@
	ln -s -f $(PWD)/$@ $(top_builddir)/include/glsim

#
# General compiler and linker options
#

AM_CPPFLAGS = -I$(top_builddir)/include
AM_CXXFLAGS=--std=c++11

#
# Distribute built sources, install headers
#

EXTRA_DIST += $(BUILT_SOURCES)

pkginclude_HEADERS = $(errors_hdr) $(logging_hdr) $(random_hdr)	\
		     $(parameters_hdr) $(environment_hdr)	\
		     $(configuration_hdr) $(simulation_hdr)	\
		     $(stochastic_hdr)

#
# binaries
#

bin_SCRIPTS = readbt.sh

#
# library
#

lib_LTLIBRARIES = libglsim.la
libglsim_la_LDFLAGS = -version-info @LIBGLSIM_VERSION@
libglsim_la_SOURCES = scontext.cc cerrors.c log.cc random.cc		\
                      parameters.cc environment.cc ctenvironment.cc	\
                      configuration.cc simulation.cc stochastic.cc
