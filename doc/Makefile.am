#
# Makefile for tex directory
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
#  glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
#                         2017, 2018, 2019
#  by Tomas S. Grigera.
# 
# glsim is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License (GPL) as published
# by the Free Software Foundation. You can use either version 3, or
# (at your option) any later version.
#
# If you use glsim to produced published work, or if you redistribute
# a modified version of glsim, or code based on glsim, please cite us
#
# - T. S. Grigera, /glsim: A general library for numerical
#   simulation/.  Computer Physics Communications *182*, 2122-2131
#   (2011).
#
# glsim is developed as part of academic work.  As such author
# attributions are not only a way to give recognition to the original
# authors, but also help to ensure continued development, since
# citations show that the work has been useful to others and thus
# help convince funding agencies that it is worth funding.  * glsim
# distribution.
#
# glsim is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.
# 
# For details see the file LICENSE in the home directory. If the file is
# missing, contact the maintainers.
#

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

doxyfile.stamp: $(top_srcdir)/prog $(top_srcdir)/lib $(top_srcdir)/lattice $(top_srcdir)/doc
	$(DOXYGEN) doxy.conf
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp

all-local: html

EXTRA_DIST = doxy.conf.in

#
# Cleaning
#

CLEANFILES += doxygen*.db
#CLEANFILES += *.dvi *.aux *.log *.toc *.nwi *.bbl *.blg glsim.ps
DISTCLEANFILES = version.tex
MAINTAINERCLEANFILES = *.tex *.def

clean-local:
	rm -rf html man latex
