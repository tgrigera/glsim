/*
 * mdobservable.cc --  Recording basic MD quantities
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
