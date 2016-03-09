/*
 * offlattice.hh -- include all offlattice headers
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, with the additional requirements of
 * attribution and nonmisrepresentation. You can use either version 3, or
 * (at your option) any later version.
 * 
 * Additional terms under GNU GPL version 3 section 7:
 * 
 * When you redistribute this software, you are required to preserve its
 * author attributions. If you distribute verbatim copies, you must not
 * alter the AUTHORS file or attributions inserted in the source files,
 * and you must not change the software's name. If you distribute a
 * modified copy, then you must give clear notice that your work is
 * different from but based on glsim. You must distribute it under a
 * different name, but include a prominent notice specifying that "(your
 * package) is based on glsim version x.x", and provide a pointer to the
 * glsim distribution.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

/** \defgroup Offlattice Off-lattice
    \brief Off-lattice library and utilities

@{

  \defgroup OfflatticeSIM Off-lattice simulation techniques

  \defgroup OfflatticeINT Description of off-lattice interactions

  \defgroup OfflatticeCONF Off-lattice configurations and trajectories

\defgroup MetricNN Metric nearest neighbours

This module includes classes to find nearest neighbours defined
metrically (i.e. those at an euclidean distance less than a certain
cutoff).

The present algorithms for nearest neighbours are to be thoght as
providing __candidates__ for NN.  For most algorithms here, pairs
returned are not guaranteed to be within cutoff; the algorithms
mereley reduce the candidates list.  In this way the support
structures need not be rebuilt at every step (which may be costly).
Candidates are guaranteed to appear only once.  For example use see
Interactions or glsim-doc.  Also check the `for_each_pair` functions
below; those will loop through all pairs and apply a given function
_only_ to those below the cutoff (i.e. they screen candidates).

@}

*/

#ifndef OFFLATTICE_HH
#define OFFLATTICE_HH

#include "ld.hh"
#include "mcobservable.hh"
#include "mdobservable.hh"
#include "lj.hh"
#include "nneighbours.hh"
#include "interactions.hh"
#include "md.hh"
#include "mc.hh"
#include "olconfiguration.hh"
#include "trajectory.hh"
#include "mdenvironment.hh"
#include "h5md.hh"

#endif /* OFFLATTICE_HH */
