/*
 * mc.cc -- Metropolis Monte Carlo for off-lattice particles with shift moves
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

#include "mc.hh"

namespace glsim {

MCParameters::MCParameters(const char *scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("MC.steps",po::value<long>()->default_value(10),"number of steps to run")
    ("MC.T",po::value<double>()->required(),"temperature")
    ("MC.DR",po::value<double>()->required(),"trial shift amplitude")
    ;
}

MCEnvironment::MCEnvironment(const char* scope) :
  SimEnvironment(scope),
  MCsteps(0),
  temperature(0),
  DR(.1),
  SE(scope),
  accepted_moves(0),
  par(scope)
{}

void MCEnvironment::init_local()
{
  SimEnvironment::init_local();
  MCsteps=par.value("MC.steps").as<long>();
  temperature=par.value("MC.T").as<double>();
  DR=par.value("MC.DR").as<double>();
  accepted_moves=0;
}

void MCEnvironment::warm_init_local()
{
  SimEnvironment::warm_init_local();
  MCsteps=par.value("MC.steps").as<long>();
  temperature=par.value("MC.T").as<double>();
  DR=par.value("MC.DR").as<double>();
  accepted_moves=0;
}

/*****************************************************************************/


MC::MC(MCEnvironment &e,OLconfiguration& c,Interactions* i) :
  Simulation(e,c),
  dran(-e.DR,e.DR),
  pran(0,1),
  env(e),
  conf(c),
  inter(i)
{
  mbeta=-1./env.temperature;
  env.total_number=conf.N;
  env.energy=0;
  update_observables();
}

void MC::step()
{
  int    i,n;
  double DE,rn[3];

  for (n=0; n<conf.N; n++) {
    rn[0]=conf.r[n][0]+dran();
    rn[1]=conf.r[n][1]+dran();
    rn[2]=conf.r[n][2]+dran();
    DE=inter->delta_energy_particle_shift(conf,n,rn);
    
    if (DE<0 || exp(mbeta*DE)>pran() ) {       /* Metropolis */
      env.accepted_moves++;
      env.energy+=DE;
      memcpy(conf.r[n],rn,sizeof(double[3]));
    }
    
  }
  inter->fold_coordinates(conf,env.DR);

  conf.step=env.steps_completed;
  conf.time=conf.step;
  env.run_completed = env.steps_in_run>=env.MCsteps;
}

void MC::update_observables()
{
  env.energy=inter->potential_energy(conf);
}

void MC::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();
  logs(info) << "   Step     Energy  Acc. rate\n";

  update_observables();
  sprintf(buff,"Initial %10.3e %10.3e\n",
	  env.steps_completed,env.energy/conf.N,
	  env.accepted_moves/(((double) env.steps_in_run)*conf.N));
  logs(info) << buff;
}

void MC::log()
{
  update_observables();
  static char buff[300];
  sprintf(buff,"%7ld %10.3e %10.3e\n",
	  env.steps_completed,env.energy/conf.N,
	  env.accepted_moves/(((double) env.steps_in_run)*conf.N));
  logs(info) << buff;
}


} /* namespace */
