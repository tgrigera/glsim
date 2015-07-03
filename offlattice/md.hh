/*
 * md.hh -- MD simulations w/different integrators
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

#ifndef MD_HH
#define MD_HH

#include "mdenvironment.hh"
#include "simulation.hh"
#include "interactions.hh"

namespace glsim {

/*****************************************************************************/

/** \class MDSimulation
    \ingroup OfflatticeSIM

Common methods for MD (does not include the integrator)

This assumes that number of particles and mass are conserved.

*/
class MDSimulation : public Simulation {
public:
  MDSimulation(MDEnvironment& e,OLconfiguration &c,Interactions *i);
  void        log();
  void        log_start_sim();

protected:
  void update_observables();

  MDEnvironment&   env;
  OLconfiguration& conf;
  Interactions     *inter;
} ;

/** \class VVerletMD
    \ingroup OfflatticeSIM

*/
class VVerletMD : public MDSimulation {
public:
  VVerletMD(MDEnvironment& e,OLconfiguration &c,Interactions *i);

  const char* name() const {return "MD with velocity Verlet integrator";}
  void step();

private:
  double           Dt,Dt2,Dtsq2;
} ;

} /* namespace */

#endif /* MD_HH */
