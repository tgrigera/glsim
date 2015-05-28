/*
 * vverlet.hh -- MD with velocity Verlet integrator
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

#ifndef VVERLETMD_HH
#define VVERLETMD_HH

#include "mdenvironment.hh"
#include "simulation.hh"
#include "interactions.hh"

namespace glsim {

/*****************************************************************************/

/** \class VVerletMD
    \ingroup Offlattice

*/
class VVerletMD : public Simulation {
public:
  VVerletMD(MDEnvironment& e,OLconfiguration &c,Interactions *i);
  ~VVerletMD() {delete[] mass;}

  const char* name() const {return "MD with velocity Verlet integrator";}
  void        step();
  void        log();
  void        log_start_sim();

  void update_observables();
  ///@{ \name Public data (for observers)
  /// Note that these quantities are only updated upon request
  double      Etot;    ///< Total energy
  double      Ekin;    ///< Kinetic energy
  double      Epot;    ///< Potential energy
  double      Ptot[3]; ///< Total momentum
  ///@}

private:
  MDEnvironment&   env;
  OLconfiguration& conf;
  Interactions     *inter;

  double           Dt,Dt2,Dtsq2;
  double           *mass;
} ;

VVerletMD::VVerletMD(MDEnvironment& e,OLconfiguration &c,Interactions *i) :
  Simulation(e,c),
  env(e),
  conf(c),
  inter(i),
  mass(0)
{
  Dt=env.time_step;
  Dt2=Dt/2;
  Dtsq2=Dt*Dt2;

  if (conf.a==0) {
    conf.a=new double[conf.N][3];
    memset(conf.a,0,conf.N*3*sizeof(double));
  }

  /// Load masses!!!!!
  mass=new double[conf.N];
  for (double *m=mass; m!=mass+conf.N; ++m)
    *m=1.;
}

} /* namespace */

#endif /* VVERLETMD_HH */
