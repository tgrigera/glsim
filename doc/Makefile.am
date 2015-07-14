#
# Makefile for tex directory
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
# glsim is copyright (C) 2009, 2010 by Tomas S. Grigera.
# 
# glsim is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License (GPL) as published by the
# Free Software Foundation, with the additional requirements of
# attribution and nonmisrepresentation. You can use either version 3, or
# (at your option) any later version.
# 
# Additional terms under GNU GPL version 3 section 7:
# 
# When you redistribute this software, you are required to preserve its
# author attributions. If you distribute verbatim copies, you must not
# alter the AUTHORS file or attributions inserted in the source files,
# and you must not change the software's name. If you distribute a
# modified copy, then you must give clear notice that your work is
# different from but based on glsim. You must distribute it under a
# different name, but include a prominent notice specifying that "(your
# package) is based on glsim version x.x", and provide a pointer to the
# glsim distribution.
#
# glsim is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.
# 
# For details see the file LICENSE in the home directory. If the file is
# missing, contact the maintainers.

#
# Suffix rules
#

# SUFFIXES = .nw .def .tex

# VPATH += : $(top_srcdir)/blib : $(top_srcdir)/lattice

#
# Targets
#

man_directory = $(top_builddir)/doc/man/man3

### tal vez dist_man_etc
##man_MANS = $(man_directory)/man_page_1.3 $(man_directory)/man_page_2.3

##$(directory)/man_page_1.3: doxyfile.stamp
##$(directory)/man_page_2.3: doxyfile.stamp

doxyfile.stamp: $(top_srcdir)/analysis $(top_srcdir)/blib $(top_srcdir)/lattice $(top_srcdir)/doc $(top_srcdir)/offlattice
	$(DOXYGEN) doxy.conf
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp

all-local: html

clean-local:
	rm -rf $(top_srcdir)/docs/man

EXTRA_DIST = doxy.conf.in

#
# Cleaning
#

CLEANFILES += doxygen*.db
#CLEANFILES += *.dvi *.aux *.log *.toc *.nwi *.bbl *.blg glsim.ps
DISTCLEANFILES = version.tex
MAINTAINERCLEANFILES = *.tex *.def

clean-local:
	rm -rf html man
