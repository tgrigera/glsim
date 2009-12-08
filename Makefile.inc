# -*- mode: makefile -*-
#
# Makefile.inc
# This include defines the directories, pattern rules and common variables
#
# version - 1.2 - 09-feb-09
# version - 1.1 - 03-dec-05
#   based on Makefile.inc v 1.2 from gsim3, 18-mar-04
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
# Copyright (C) 2009 by Tomas S. Grigera.
#
# glsim is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# glsim is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.
#
# For details, see the file LICENSE in the base directory. If the
# file is missing, look for the GNU General Public License v3 in
# <http://www.gnu.org/licenses/>.
#

#
# Directories
#

INCDIRS=-I$(glsimd)/base -I/usr/include/netcdf

#
# Compiler flags (target- and compiler-specific)
#
CPPFLAGS=$(INCDIRS)
LIBDIRS=$(patsubst %,-L%,$(libd))
INSTALLD=$(bind)

include $(glsimd)/Makefile.options

#
# Libraries
#

FORTRANLIBS=-lgfortran
NETCDFLIBS=-lnetcdf_c++ -lnetcdf
BOOSTLIBS=-lboost_program_options
GSL=-lgsl -lgslcblas

COMMON_3DPARTY=$(NETCDFLIBS) $(BOOSTLIBS) $(GSL)
LIBRARIES=$(COMMON_3DPARTY) $(SYSLIBS)

#
# Pattern rules we use (I must add the target directory)
#

TANGLE=notangle
WEAVE=noweave

%.c:
	$(TANGLE) -L -R$@ $< | cpif $@

%.cc:
	$(TANGLE) -L -R$@ $< | cpif $@

%.cdl:
	$(TANGLE) -L -R$@ $< | cpif $@

%.f:
	$(TANGLE) -L -R$@ $< | cpif $@

%.h:
	$(TANGLE) -L -R$@ $< | cpif $@

%.hh:
	$(TANGLE) -L -R$@ $< | cpif $@

%.sh:
	$(TANGLE) -L -R$@ $< | cpif $@

$(TDIR)/%.awk:
	$(TANGLE) -L -R$(@F) $< >$(@F)
	chmod +x $(@F)
	ln -s -f $$PWD/$(@F) $@

$(TDIR)/%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

$(TDIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ -c $<

$(TDIR)/%.o: %.f
	$(FC) $(FFLAGS) -o $@ -c $<

%.dd: %.c
	gcc $(CPPFLAGS) -MM -MT '$$(TDIR)/$*.o' $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%.dd: %.cc
	gcc $(CPPFLAGS) -MM -MT '$$(TDIR)/$*.o' $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%.dd: %.f
	echo '$$(TDIR)/$*.o $*.dd: $<' > $@

#
# In the .d's I want to include the dependencies as well as the object list
#

%.d: %.nw
	@echo "$@: $^" >$@ ; \
	SRCS=$$(noroots $^ | sed 's/<<\(.*\)>>/\1/' ) ; \
	echo $$SRCS | tr \\n ' ' | sed 's/\(.*\)/\1: $^\n/' >>$@ ; \
	echo "sources+=$$SRCS" | tr \\n ' ' >>$@ ; \
	echo " " >>$@ ; \
	echo "$$SRCS" | tr ' ' \\n | egrep '*.\.[cc|c|f]' | \
	  sed 's|\(.*\.\)\(.*\)|\1dd: \1\2 \nobjects+=$$(TDIR)/\1o\ninclude \1dd|' >>$@
#	echo "$$SRCS" | tr ' ' \\n | egrep '*.\.[cc|c|f]' | \
#	  sed 's|\(.*\.\)\(.*\)|\1dd: \1\2 $$(sources)\nobjects+=$$(TDIR)/\1o\ninclude \1dd|' >>$@

$(TDIR)/%::
	$(LINKER_CALLER) $(LDFLAGS) -o $@ $+ $(LIBRARIES)

#
# To read list of source files automatically
#

websources := $(wildcard *.nw)
depends = $(patsubst %.nw,%.d,$(websources))

.SUFFIXES:

# Search paths

include $(glsimd)/Makefile.libpaths
VPATH=$(TDIR)
GPATH=$(VPATH)

execs=$(patsubst %,$(TDIR)/%,$(exec))

.PHONY: all

# I need all twice: the first time to make it the first target, then again
# to specify the dependencies after $(depends) has been loaded

all:

include $(depends)

all: $(sources) $(objects) $(execs)
