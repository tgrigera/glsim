/*
 * ld.cc -- Langevin dynamics simulations
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

#include "ld.hh"

namespace glsim {

/*****************************************************************************
 *
 * Parameters and environment
 *
 */
  
LDParameters::LDParameters(const char *scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("LD.T",po::value<double>()->required(),"temperature")
    ("LD.eta",po::value<double>()->required(),"friction coefficient")
    ;
}

LDEnvironment::LDEnvironment(const char* scope) :
  MDEnvironment(scope),
  eta(0),
  temperature(0),
  SE(scope),
  par(scope)
{}

void LDEnvironment::init_local()
{
  MDEnvironment::init_local();
  temperature=par.value("LD.T").as<double>();
  eta=par.value("LD.eta").as<double>();
}

void LDEnvironment::warm_init_local()
{
  MDEnvironment::warm_init_local();
  temperature=par.value("LD.T").as<double>();
  eta=par.value("LD.eta").as<double>();
}

template <typename Archive>
inline void LDEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw glsim::Environment_wrong_version("MDEnvironment",version,class_version);
  ar & boost::serialization::base_object<MDEnvironment>(*this);
  ar & temperature & eta;
}

/*****************************************************************************
 *
 * Simulation
 *
 */

LDSimulation::LDSimulation(LDEnvironment& e,OLconfiguration &c,Interactions *i) :
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

  // Constants for Langevin integration
  Dt=env.time_step;
  for (int i=0; i<conf.N; ++i)
    if (c0.count(conf.type[i])==0) {
      double mass=inter->mass(conf.type[i]);
      double xi=env.eta/mass;
      double xidt=xi*Dt;
      double c0l,c1,c2,sx,sv,rho;
      double exi=exp(-xidt);
      if (xidt<1e-3) {
	c0l=1 - xidt + xidt*xidt/2 - xidt*xidt*xidt/6;
	c1=1 - xidt/2 + xidt*xidt/6 - xidt*xidt*xidt/24;
	c2=0.5 - xidt/6 + xidt*xidt/24;
      } else {
	c0l=exi;
	c1=(1-c0l)/xidt;
        c2=(1-c1)/xidt;
      }
      sv=(env.temperature/inter->mass(conf.type[i]))*(1-exi*exi);
      sv=sqrt(sv);
      if (env.eta<1e-3 || xi<1e-4) {
	sx=env.temperature*Dt*Dt*Dt*env.eta*(2./3.-0.5*xidt)/(mass*mass);
	sx=sqrt(sx);
	rho=sqrt(3.)*(0.5-xidt/16.-(17./1280.)*xidt*xidt
		     +(17./6144)*xidt*xidt*xidt);
      } else {
	sx=(env.temperature/env.eta)*(2*Dt-(3-4*exi+exi*exi)/xi);
	sx=sqrt(sx);
	double sxv=(env.temperature/env.eta)*(1-exi)*(1-exi);
	rho=sxv/(sx*sv);
      }
      	
      c0[conf.type[i]]=c0l;
      c1dt[conf.type[i]]=Dt*c1;
      c2dtsq[conf.type[i]]=Dt*Dt*c2;
      c1mc2dt[conf.type[i]]=(c1-c2)*Dt;
      c2dt[conf.type[i]]=c2*Dt;
      noise[conf.type[i]]=new BivariateGaussian_distribution(sx,sv,rho);
    }
}

LDSimulation::~LDSimulation()
{
  for (auto &p : noise)
    delete p.second;
}

void LDSimulation::update_observables()
{
  env.Ekin=inter->kinetic_energy_and_momentum(conf,env.Ptot);
  env.Etot=env.Ekin+env.Epot;
}

void LDSimulation::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();
  logs(info) << "    Step       Time         Px         Py         Pz       Epot       Ekin       Etot         kT\n";

  sprintf(buff," Initial            %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.Ptot[0],env.Ptot[1],env.Ptot[2],env.Epot/conf.N,env.Ekin/conf.N,
	  env.Etot/conf.N,(2./3.)*(env.Ekin/conf.N));
  logs(info) << buff;
}

void LDSimulation::log()
{
  update_observables();
  static char buff[300];
  sprintf(buff,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	  env.Ptot[0],env.Ptot[1],env.Ptot[2],
	  env.Epot/conf.N,env.Ekin/conf.N,env.Etot/conf.N,(2./3)*(env.Ekin/conf.N));
  logs(info) << buff;
}

void LDSimulation::step()
{
  double maxdisp=0;
  double drx,dry,drz,xir,xiv;

  for (int i=0; i<conf.N; ++i) {
    short t=conf.type[i];
    (*noise[t])(xir,xiv);
    conf.r[i][0] += ( drx = c1dt[t]*conf.v[i][0] + c2dtsq[t]*conf.a[i][0] + xir );
    conf.v[i][0] = c0[t]*conf.v[i][0] + c1mc2dt[t]*conf.a[i][0] + xiv;
    (*noise[t])(xir,xiv);
    conf.r[i][1] += ( dry = c1dt[t]*conf.v[i][1] + c2dtsq[t]*conf.a[i][1] + xir );
    conf.v[i][1] = c0[t]*conf.v[i][1] + c1mc2dt[t]*conf.a[i][1] + xiv;
    (*noise[t])(xir,xiv);
    conf.r[i][2] += ( drz = c1dt[t]*conf.v[i][2] + c2dtsq[t]*conf.a[i][2] + xir );
    conf.v[i][2] = c0[t]*conf.v[i][2] + c1mc2dt[t]*conf.a[i][2] + xiv;

    double drsq=drx*drx + dry*dry + drz*drz;
    if (drsq>maxdisp) maxdisp=drsq;
  }
  maxdisp=sqrt(maxdisp);  // This is the modulus of the largest
			  // displacement, to be used by
			  // fold_coordinates to rebuild neighbour
			  // list when needed
  env.Epot=inter->acceleration_and_potential_energy(conf);
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0] += c2dt[conf.type[i]]*conf.a[i][0];
    conf.v[i][1] += c2dt[conf.type[i]]*conf.a[i][1];
    conf.v[i][2] += c2dt[conf.type[i]]*conf.a[i][2];
  }
  inter->fold_coordinates(conf,maxdisp);

  env.time_completed+=Dt;
  env.time_in_run+=Dt;
  env.run_completed = env.steps_in_run>=env.MDsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
}

} /* namespace */
