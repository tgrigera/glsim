/*
 * vverlet.cc -- MD with velocity Verlet integrator
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

#include "vverletmd.hh"

namespace glsim {

VVerletMD::VVerletMD(MDEnvironment& e,OLconfiguration &c,Interactions *i) :
  Simulation(e,c),
  env(e),
  conf(c),
  inter(i)
{
  Dt=env.time_step;
  Dt2=Dt/2;
  Dtsq2=Dt*Dt2;

  if (conf.a==0) {
    conf.a=new double[conf.N][3];
    memset(conf.a,0,conf.N*3*sizeof(double));
  }
  if (conf.v==0) {
    conf.v=new double[conf.N][3];
    // maxwell!!!
  }
  // Substract Vcm
  total_mass=0;
  for (int i=0; i<conf.N; ++i) {
    total_mass+=inter->mass(conf.type[i]);
    Ptot[0]+=conf.v[i][0]*inter->mass(conf.type[i]);
    Ptot[1]+=conf.v[i][1]*inter->mass(conf.type[i]);
    Ptot[2]+=conf.v[i][2]*inter->mass(conf.type[i]);
  }
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0]-=Ptot[0]/total_mass;
    conf.v[i][1]-=Ptot[1]/total_mass;
    conf.v[i][2]-=Ptot[2]/total_mass;
  }

  Epot=inter->potential_energy(conf);
  update_observables();
}

void VVerletMD::step()
{
  for (int i=0; i<conf.N; ++i) {
    conf.r[i][0] += Dt*conf.v[i][0] + Dtsq2*conf.a[i][0];
    conf.r[i][1] += Dt*conf.v[i][1] + Dtsq2*conf.a[i][1];
    conf.r[i][2] += Dt*conf.v[i][2] + Dtsq2*conf.a[i][2];
    conf.v[i][0] += Dt2*conf.a[i][0];
    conf.v[i][1] += Dt2*conf.a[i][1];
    conf.v[i][2] += Dt2*conf.a[i][2];
  }
  Epot=inter->acceleration_and_potential_energy(conf);
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0] += Dt2*conf.a[i][0];
    conf.v[i][1] += Dt2*conf.a[i][1];
    conf.v[i][2] += Dt2*conf.a[i][2];
  }
  inter->fold_coordinates(conf);

  env.time_completed+=Dt;
  env.time_in_run+=Dt;
  env.run_completed = env.steps_in_run>=env.MDsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
}

void VVerletMD::update_observables()
{
  Ekin=inter->kinetic_energy_and_momentum(conf,Ptot);
  Etot=Ekin+Epot;
}

void VVerletMD::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();
  logs(info) << "   Step       Time         Px         Py         Pz       Epot       Ekin       Etot         kT\n";

  sprintf(buff,"Initial            %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  Ptot[0],Ptot[1],Ptot[2],Epot/conf.N,Ekin/conf.N,Etot/conf.N,
	  (2./3.)*(Ekin/conf.N));
  logs(info) << buff;
}

void VVerletMD::log()
{
  update_observables();
  static char buff[300];
  sprintf(buff,"%7ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,Ptot[0],Ptot[1],Ptot[2],
	  Epot/conf.N,Ekin/conf.N,Etot/conf.N,(2./3)*(Ekin/conf.N));
  logs(info) << buff;
}

} /* namespace */
