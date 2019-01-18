/*
 * analysis.hh -- Doxygen page for analysis directory
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use glsim to produced published work, or if you redistribute a
 * modified version of glsim, or code based on glsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * glsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding.  * glsim
 * distribution.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

/** @{
    \defgroup Analysis Analysis
    \brief Library and programs to analyse simulation data

    @{    

    \defgroup Structure Structure
     \brief    Utilities and classes to compute structural properties of off-lattice systems


    @}

    @{
    \defgroup TCorr Time correlations

     The programs and routines in this module compute a discretized
     two-point time correlation function. Assuming time-traslation
     invariance (i.e.\ independence of the time origin), the
     correlation can be computed as

     \f[  C(t) = \sum_{t_0} A(t+t_0)  A^*(t_0) \mbox{NORMALIZATION MISSING}.\f]

     The connected correlation is defined as

     \f[$ C_c(t) = C(t) - \langle A     \rangle^2. \f]

    @}

    @}
*/

