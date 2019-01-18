/*
 * mdobservable.cc --  Recording basic MD quantities
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

#include "mdobservable.hh"

namespace glsim {

MDObservable_parameters::MDObservable_parameters(const char* scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("MD.obs_interval",po::value<int>()->default_value(0),
     "Interval for standard observation, 0=skip")
    ("MD.obs_file_prefix",po::value<std::string>(),"Observation file prefix")
    ;
}

void MDObservable::interval_and_file()
{
  obs_interval=par.value("MD.obs_interval").as<int>();
  obs_file_prefix=par.value("MD.obs_file_prefix").as<std::string>();
}

void MDObservable::write_header()
{
  fprintf(of,"    Step       Time         Px         Py         Pz       Epot       Ekin       Etot         kT\n");
}

void MDObservable::observe()
{
  env.Ekin=env.inter->kinetic_energy_and_momentum(conf,env.Ptot);
  env.Etot=env.Ekin+env.Epot;
  int N=env.total_number;
  fprintf(of,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	  env.Ptot[0],env.Ptot[1],env.Ptot[2],
	  env.Epot/N,env.Ekin/N,env.Etot/N,(2./3)*(env.Ekin/N));
}

} /* namespace */
