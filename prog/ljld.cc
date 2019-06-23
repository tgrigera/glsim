/*
 * ljld.hh -- Langevin Dynamics of Lennard-Jones particles
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
 * help convince funding agencies that it is worth funding.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

/** \file ljld.cc
    \ingroup Offlattice
    \brief Langevin dynamics of repulsive Lennard-Jones particles
*/

#include "olconfiguration.hh"
#include "mdenvironment.hh"
#include "mdobservable.hh"
#include "trajectory.hh"
#include "lj.hh"
#include "ld.hh"

void wmain(int argc, char *argv[])
{
  glsim::LDEnvironment  env;
  glsim::OLconfiguration conf;
  glsim::RepulsiveLennardJones LJ(env.scope());
  glsim::MDObservable obs(env,conf);
  glsim::SimulationCL CL("GS_ljmd","(C) 2015 Tomas S. Grigera",env.scope());
  glsim::Trajectory traj(env,conf,glsim::OLconfig_file::options().r_frame());
  
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  glsim::Interactions_isotropic_pairwise<glsim::RepulsiveLennardJones>
    inter(LJ,conf);
  glsim::LDSimulation sim(env,conf,&inter);
  obs.observe_first();   // This is optional
  traj.observe_first();  // This is mandatory!
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
