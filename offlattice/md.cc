/*
 * md.cc -- MD simulations
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

#include "md.hh"

namespace glsim {

MDSimulation::MDSimulation(MDEnvironment& e,OLconfiguration &c,Interactions *i) :
  Simulation(e,c),
  env(e),
  conf(c),
  inter(i)
{
  if (conf.a==0) {
    conf.a=new double[conf.N][3];
    memset(conf.a,0,conf.N*3*sizeof(double));
  }
  if (conf.v==0) {
    conf.v=new double[conf.N][3];
    memset(conf.v,0,conf.N*sizeof(double));
    // maxwell!!!
  }

  env.total_number=conf.N;
  env.total_mass=0;
  env.inter=inter;
  for (int i=0; i<conf.N; ++i) {
    env.total_mass+=inter->mass(conf.type[i]);
    env.Ptot[0]+=conf.v[i][0]*inter->mass(conf.type[i]);
    env.Ptot[1]+=conf.v[i][1]*inter->mass(conf.type[i]);
    env.Ptot[2]+=conf.v[i][2]*inter->mass(conf.type[i]);
  }
  // Substract Vcm (if P is conserved)
  if (inter->conserve_P())
    for (int i=0; i<conf.N; ++i) {
      conf.v[i][0]-=env.Ptot[0]/env.total_mass;
      conf.v[i][1]-=env.Ptot[1]/env.total_mass;
      conf.v[i][2]-=env.Ptot[2]/env.total_mass;
    }

  env.Epot=inter->potential_energy(conf);
  conf.step=env.steps_completed;
  conf.time=env.time_completed;
  update_observables();
}

void MDSimulation::update_observables()
{
  env.Ekin=inter->kinetic_energy_and_momentum(conf,env.Ptot);
  env.Etot=env.Ekin+env.Epot;
}

void MDSimulation::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();
  logs(info) << "    Step       Time         Px         Py         Pz       Epot       Ekin       Etot         kT\n";

  sprintf(buff," Initial            %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.Ptot[0],env.Ptot[1],env.Ptot[2],env.Epot/conf.N,env.Ekin/conf.N,
	  env.Etot/conf.N,(2./3.)*(env.Ekin/conf.N));
  logs(info) << buff;
}

void MDSimulation::log()
{
  update_observables();
  static char buff[300];
  sprintf(buff,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	  env.Ptot[0],env.Ptot[1],env.Ptot[2],
	  env.Epot/conf.N,env.Ekin/conf.N,env.Etot/conf.N,(2./3)*(env.Ekin/conf.N));
  logs(info) << buff;
}

/*****************************************************************************
 *
 * VVerletMD
 *
 */

VVerletMD::VVerletMD(MDEnvironment& e,OLconfiguration &c,Interactions *i) :
  MDSimulation(e,c,i)
{
  Dt=env.time_step;
  Dt2=Dt/2;
  Dtsq2=Dt*Dt2;
}

void VVerletMD::step()
{
  double maxdisp=0;
  double drx,dry,drz;

  for (int i=0; i<conf.N; ++i) {
    conf.r[i][0] += ( drx = Dt*conf.v[i][0] + Dtsq2*conf.a[i][0] );
    conf.r[i][1] += ( dry = Dt*conf.v[i][1] + Dtsq2*conf.a[i][1] );
    conf.r[i][2] += ( drz = Dt*conf.v[i][2] + Dtsq2*conf.a[i][2] );
    conf.v[i][0] += Dt2*conf.a[i][0];
    conf.v[i][1] += Dt2*conf.a[i][1];
    conf.v[i][2] += Dt2*conf.a[i][2];

    double drsq=drx*drx + dry*dry + drz*drz;
    if (drsq>maxdisp) maxdisp=drsq;
  }
  maxdisp=sqrt(maxdisp);  // This is the modulus of the largest
			  // displacement, to be used by
			  // fold_coordinates to rebuild neighbour
			  // list when needed
  env.Epot=inter->acceleration_and_potential_energy(conf);
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0] += Dt2*conf.a[i][0];
    conf.v[i][1] += Dt2*conf.a[i][1];
    conf.v[i][2] += Dt2*conf.a[i][2];
  }
  inter->fold_coordinates(conf,maxdisp);

  env.time_completed+=Dt;
  env.time_in_run+=Dt;
  env.run_completed = env.steps_in_run>=env.MDsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
}

} /* namespace */
