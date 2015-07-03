#
# Makefile.am for offlattice directory
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
# glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015
# by Tomas S. Grigera.
#
# See file LICENSE in base directory
#

#
# General compiler and linker options
#

include $(top_srcdir)/Makefile-common
AM_CPPFLAGS = -I$(top_srcdir)/blib

#
# Headers
#

pkginclude_HEADERS = h5md.hh olconfiguration.hh interactions.hh ld.hh	\
		     lj.hh mc.hh mdenvironment.hh mdobservable.hh	\
		     md.hh nneighbours.hh trajectory.hh

#
# Offlattice library
#

lib_LTLIBRARIES = libglsim_ol.la
libglsim_ol_la_LDFLAGS = -version-info @LIBGLSIM_VERSION@
libglsim_ol_la_SOURCES = h5md.cc olconfiguration.cc interactions.cc	\
			 ld.cc mc.cc mcobservable.cc md.cc		\
			 mdenvironment.cc mdobservable.cc		\
			 nneighbours.cc trajectory.cc

#
# Binaries
#

LDADD=$(top_builddir)/blib/libglsim.la libglsim_ol.la

bin_PROGRAMS = GS_olconf_cat GS_olconf_create GS_ljld GS_ljmc GS_ljmd GS_traj2h5md

GS_olconf_cat_SOURCES = olconf_cat.cc
GS_olconf_create_SOURCES = olconf_create.cc
GS_ljld_SOURCES = ljld.cc
GS_ljmc_SOURCES = ljmc.cc
GS_ljmd_SOURCES = ljmd.cc
GS_traj2h5md_SOURCES = traj2h5md.cc

#
# Test programs
#

check_PROGRAMS = pseudo-md

pseudo_md_SOURCES = pseudo-md.cc


