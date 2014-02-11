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

SUFFIXES = .nw .def .tex

VPATH += : $(top_srcdir)/blib

.nw.def:
	@DEF@ $< >$@

.nw.tex:
	@WEAVE@ -delay -indexfrom all.def -filter 'sed "/^@use / s/_/\\\\_/g; /^@defn / s/_/\\\\_/g;"' $< >$@

#
# Targets
#

BUILT_SOURCES = glsim.dvi

glsim.dvi: $(texfiles)
	export TEXINPUTS=$(TEXINPUTS):$(top_srcdir)/doc: ; \
	export BIBINPUTS=$(top_srcdir)/doc ; \
        @LATEX@ glsim.tex ; \
	@BIBTEX@ glsim ; \
	@INDEX@ glsim.tex ; \
	sed 's/_/\\_/g;' glsim.nwi >glsim.nwi.$$$$ ; mv glsim.nwi.$$$$ glsim.nwi ;\
        @LATEX@ glsim.tex ; \
        @LATEX@ glsim.tex

texfiles = glsim.tex $(blib_texfiles)

deffiles = glsim.def $(blib_deffiles)

blib_texfiles = blib.tex errors.tex logging.tex random.tex parameters.tex

blib_deffiles = blib.def errors.def logging.def random.def parameters.def

$(texfiles): all.def

all.def: $(deffiles)
	sort -u $^ >$@

#
# Installation
#

#
# Distribution
#

EXTRA_DIST = glsim.nw glsim.bib gpl.tex version.tex.in $(texfiles) \
             all.def $(deffiles)

#
# Cleaning
#

CLEANFILES = *.dvi *.aux *.log *.toc *.nwi *.bbl *.blg glsim.ps
DISTCLEANFILES = version.tex
MAINTAINERCLEANFILES = *.tex *.def
