/*
 * olglsim.hh -- include all headers
 *
 * This file is part of olglsim, a numerical simulation class library
 * and helper programs.
 *
 * olglsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * olglsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use olglsim to produced published work, or if you redistribute a
 * modified version of olglsim, or code based on olglsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * olglsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding.
 *
 * olglsim is distributed in the hope that it will be useful, but WITHOUT
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
