/*
 * mcobservable.cc --  Recording basic quantities along a MC run
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

#include "mcobservable.hh"

namespace glsim {

MCObservable_parameters::MCObservable_parameters(const char* scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("MC.obs_interval",po::value<int>()->default_value(0),
     "Interval for standard observation, 0=skip")
    ("MC.obs_file_prefix",po::value<std::string>(),"Observation file prefix")
    ;
}

void MCObservable::interval_and_file()
{
  obs_interval=par.value("MC.obs_interval").as<int>();
  obs_file_prefix=par.value("MC.obs_file_prefix").as<std::string>();
}

void MCObservable::write_header()
{
  fprintf(of,"    Step     Energy   Acc.Rate\n");
}

void MCObservable::observe()
{
  int N=env.total_number;
  fprintf(of,"%8ld %10.3e %10.3e\n",
	  env.steps_completed,
	  env.energy/N,(double) env.accepted_moves/(N*env.steps_in_run));
}

} /* namespace */
