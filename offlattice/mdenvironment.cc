/*
 * mdenvironment.cc -- Environment for Molecular Dynamics -- definitions
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

#include "mdenvironment.hh"

namespace glsim {

/*****************************************************************************
 *
 * MDParameters
 *
 */

MDParameters::MDParameters(const char *scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("MD.steps",po::value<long>()->default_value(10),"number of steps to run")
    ("MD.time_step",po::value<double>()->required(),"time step for integrator")
    ;
}

/*****************************************************************************
 *
 * MDEnvironmet
 *
 */

MDEnvironment::MDEnvironment(const char* scope) :
  SimEnvironment(scope),
  MDsteps(0),
  time_step(0.),
  par(scope)
{}

void MDEnvironment::init_local()
{
  SimEnvironment::init_local();
  MDsteps=par.value("MD.steps").as<long>();
  time_step=par.value("MD.time_step").as<double>();
}

void MDEnvironment::warm_init_local()
{
  SimEnvironment::warm_init_local();
  MDsteps=par.value("MD.steps").as<long>();
  time_step=par.value("MD.time_step").as<double>();
}

template <typename Archive>
inline void MDEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw glsim::Environment_wrong_version("MDEnvironment",version,class_version);
  ar & boost::serialization::base_object<SimEnvironment>(*this);
  ar & MDsteps;
  ar & time_step;
}

} /* namespace */
