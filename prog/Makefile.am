#
# Makefile.am for analysis directory
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
# General compiler and linker options
#

include $(top_srcdir)/Makefile-common
AM_CPPFLAGS = -I$(top_srcdir)/lib/blib -I$(top_srcdir)/lib/offlattice -I$(top_srcdir)/lib/analysis

#
# binaries
#

bin_SCRIPTS = gs_columnave.awk gs_fileave.awk readbt.sh

EXTRA_DIST = readbt.sh gs_columnave.awk gs_fileave.awk

LDADD=$(top_builddir)/lib/libglsim.la

bin_PROGRAMS = gs_avesd gs_gav gs_histogram gs_tcorr gs_tau

gs_avesd_SOURCES = avesd.cc
gs_gav_SOURCES = gav.cc
gs_histogram_SOURCES = histocl.cc
gs_tcorr_SOURCES = tcorr.cc
gs_tau_SOURCES = tau.cc

#
# Offlattice progs
#

if HAS_HDF

LDADD+=$(top_builddir)/lib/libglsim_ol.la

bin_PROGRAMS += gs_1stneighbours gs_gr gs_Fk gs_Sk gs_olconf_cat	\
		gs_olconf_create gs_olconf_dump gs_ljenergy gs_ljld	\
		gs_ljmc gs_ljmd gs_neighbour_time gs_rljmc gs_rljmd	\
		gs_traj2h5md

gs_1stneighbours_SOURCES = gs_1stneighbours.cc
gs_1stneighbours_CXXFLAGS=$(OPENMP_CXXFLAGS) $(AM_CXXFLAGS)
gs_Fk_SOURCES = Fk.cc
gs_gr_SOURCES = gs_gr.cc
gs_Sk_SOURCES = gs_Sk.cc
gs_olconf_cat_SOURCES = olconf_cat.cc
gs_olconf_create_SOURCES = olconf_create.cc
gs_olconf_dump_SOURCES = olconf_dump.cc
gs_ljenergy_SOURCES = ljenergy.cc energy.cc
gs_ljld_SOURCES = ljld.cc
gs_ljmc_SOURCES = ljmc.cc
gs_ljmd_SOURCES = ljmd.cc
gs_neighbour_time_SOURCES = neighbour_time.cc
gs_neighbour_time_CXXFLAGS=$(OPENMP_CXXFLAGS) $(AM_CXXFLAGS)
gs_rljmc_SOURCES = rljmc.cc
gs_rljmd_SOURCES = rljmd.cc
gs_traj2h5md_SOURCES = traj2h5md.cc

endif
