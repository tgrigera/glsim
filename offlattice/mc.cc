/*
 * mc.cc -- Metropolis Monte Carlo for off-lattice particles with shift moves
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

template <typename Archive>
inline void MCEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw glsim::Environment_wrong_version("MDEnvironment",version,class_version);
  ar & boost::serialization::base_object<SimEnvironment>(*this);
  ar & MCsteps;
  ar & temperature & DR;
  ar & accepted_moves;
}


/*****************************************************************************/


MC::MC(MCEnvironment &e,OLconfiguration& c,Interactions* i) :
  Simulation(e,c),
  energy(0),
  dran(-e.DR,e.DR),
  pran(0,1),
  env(e),
  conf(c),
  inter(i)
{
  mbeta=-1./env.temperature;
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
      energy+=DE;
      memcpy(conf.r[n],rn,sizeof(double[3]));
    }
    
  }
  inter->fold_coordinates(conf);

  conf.step=env.steps_completed;
  conf.time=conf.step;
  env.run_completed = env.steps_in_run>=env.MCsteps;
}

void MC::update_observables()
{
  energy=inter->potential_energy(conf);
}

void MC::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();
  logs(info) << "   Step     Energy  Acc. rate\n";

  update_observables();
  sprintf(buff,"Initial %10.3e %10.3e\n",
	  env.steps_completed,energy/conf.N,
	  env.accepted_moves/(((double) env.steps_in_run)*conf.N));
  logs(info) << buff;
}

void MC::log()
{
  update_observables();
  static char buff[300];
  sprintf(buff,"%7ld %10.3e %10.3e\n",
	  env.steps_completed,energy/conf.N,
	  env.accepted_moves/(((double) env.steps_in_run)*conf.N));
  logs(info) << buff;
}


} /* namespace */
