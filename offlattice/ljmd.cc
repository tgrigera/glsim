/*
 * ljmd.hh -- Molecular Dynamics of Lennard-Jones particles
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

/** \file ljmd.cc
    \ingroup Offlattice
    \brief Molecular dynamics of Lennard-Jones particles
*/

#include "olconfiguration.hh"
#include "mdenvironment.hh"
#include "mdobservable.hh"
#include "trajectory.hh"
#include "lj.hh"
#include "md.hh"

void wmain(int argc, char *argv[])
{
  glsim::MDEnvironment  env;
  glsim::OLconfiguration conf;
  glsim::LennardJones LJ(env.scope());
  glsim::MDObservable obs(env,conf);
  glsim::SimulationCL CL("GS_ljmd","(C) 2015 Tomas S. Grigera",env.scope());
  glsim::Trajectory traj(env,conf,glsim::OLconfig_file::options().r_frame());
  
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  glsim::Interactions_isotropic_pairwise<glsim::LennardJones>
    inter(LJ,conf);
  glsim::VVerletMD sim(env,conf,&inter);
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
