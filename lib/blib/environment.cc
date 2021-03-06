/*
 * environment.cc -- definitions for environment classes
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use glsim to produced published work, or if you redistribute a
 * modified version of glsim, or code based on glsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * glsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding.  * glsim
 * distribution.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

#include <fstream>
#include <float.h>

#include "environment.hh"

namespace glsim {

/******************************************************************************
 *
 * Environment
 *
 */

// The default constructor below adds the object to the list of
// environments belonging to the specified scope.  This is managed with a
// [[std::map]] object.  The destructor clears it from the list.
Environment::scopemap_t  Environment::scopes;

Environment::Environment(const char* scope) : 
  initialized(deflt),
  scope_name(scope),
  ini_infix("+++"),
  fin_infix("xxx"),
  extension(".dat"),
  par(scope)
{
#ifdef DEBUG
  logs(debug) << "(DD) (Environment) Pushing to scope " << scope << " ("
	      << this << ")\n";
#endif
  scopes[scope_name].push_back(this);
}

Environment::~Environment()
{
  scopes[scope_name].remove(this);
}

// Initialization: in this case there is no difference between full
// and warm, so warm_init() just calls init(). Note however that we
// must be careful because init() and warm_init() are virtual.  At
// this level, it is correct to implement warm_init() by calling
// init(), but if we must avoid the virtual call mechanism (by
// qualifying it with the class name), because otherwise way we would
// in effect be invoking the init() of the __derived__ class, which is
// not what we want.

// Note also that ini_infix and fin_infix are expected from a
// parameters object to be declared elsewhere.  Normally this would be
// a StandardCL with the same scope, but you may provide otherwise.
// If you do not wish to use the filename creation functions just
// ignore them, it is not an error if the infixes are not declared.

void Environment::warm_init_local()
{
  initialized=warm;
  Environment::init_common();
}

void Environment::init_local()
{
  initialized=cold;
  Environment::init_common();
}

void Environment::init_common()
{
  if (par.count("initial_infix")>0)
      ini_infix=par.value("initial_infix").as<std::string>();
  if (par.count("final_infix")>0)
    fin_infix=par.value("final_infix").as<std::string>();
}

// Finally, these methods can be used to build an initial and final
// filenames.  If the first argument contains the string [["[AUTO]"]],
// then it will be replaced with a name of the form
// prefix+infix+extension.  The extension is [[.dat]] by default, but you
// can change it initializing [[extension]] to something else in your
// derived [[init()]] method.

void Environment::initial_filename(std::string& name,
					const std::string& prefix)
{
  if (name!="[AUTO]") return;
  name=prefix+ini_infix+extension;
}

void Environment::final_filename(std::string& name,
				      const std::string& prefix)
{
  if (name!="[AUTO]") return;
  name=prefix+fin_infix+extension;
}

/******************************************************************************
 *
 * BaseEnv
 *
 */

BaseEnv_par::BaseEnv_par(const char* scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("initial_configuration_file",po::value<std::string>()->default_value("[AUTO]"),
     "file to read inital configuration, if [AUTO], generated from prefix and infix")
    ("final_configuration_file",po::value<std::string>()->default_value("[AUTO]"),
     "file to write final configuration, if [AUTO], generated from prefix and infix")
    ("initial_environment_file",po::value<std::string>()->default_value("[AUTO]"),
     "file to read inital environment, if [AUTO], generated from prefix and infix")
    ("final_environment_file",po::value<std::string>()->default_value("[AUTO]"),
     "file to write final environment, if [AUTO], generated from prefix and infix")
    ("configuration_file_prefix",po::value<std::string>()->default_value("conf"),
     "prefix to generate initial and final configuration filenames")
    ("environment_file_prefix",po::value<std::string>()->default_value("env"),
     "prefix to generate initial and final environment filenames")
    ;
}

//  Now the (object-level) initialization functions.  The local
// initizializtion is public because partial initialization is useful
// when checking for a partially completed simulation.  Note that the
// requirement to call the base class init functions is fulfilled.

void BaseEnvironment::init_base()
{
  BaseEnvironment::init_local();
}

void BaseEnvironment::warm_init_local()
{
  Environment::warm_init_local();
  init_local_common();
}

void BaseEnvironment::init_local()
{
  Environment::init_local();
  init_local_common();
}

void BaseEnvironment::init_local_common()
{
  configuration_file_ini=par.value("initial_configuration_file").as<std::string>();
  configuration_file_fin=par.value("final_configuration_file").as<std::string>();
  environment_file_ini=par.value("initial_environment_file").as<std::string>();
  environment_file_fin=par.value("final_environment_file").as<std::string>();
  configuration_file_prefix=par.value("configuration_file_prefix").as<std::string>();
  environment_file_prefix=par.value("environment_file_prefix").as<std::string>();

  initial_filename(configuration_file_ini,configuration_file_prefix);
  final_filename(configuration_file_fin,configuration_file_prefix);
  initial_filename(environment_file_ini,environment_file_prefix);
  final_filename(environment_file_fin,environment_file_prefix);
}

// @ \subsection{Scope-level functions}

// Parsing the scope is actually implemented in [[Parameters]].
// Initialization of the scope simply involves going through the list of
// environments and calling the corresponding methods.
void BaseEnvironment::init()
{
  scope_t& scope=scopes[scope_name];

  for (Environment* e : scope) e->init_local();
}

void BaseEnvironment::warm_init()
{
  scope_t& scope=scopes[scope_name];

  for (Environment* e : scope) e->warm_init_local();
}

void BaseEnvironment::parse(const char* file)
{
  par.parse(file);
}

void BaseEnvironment::step()
{
  scope_t& scope=scopes[scope_name];

  for (Environment* e : scope) e->step_local();
}

// Loading and saving call the serialization methods after opening a
// binary archive.  Note that serialization here is done with
// \emph{objects} (actually, references) and not pointers.  Thus the
// environments must be created \emph{before} loading.  This is usually
// the case in the simulation use, but you may want some environments to
// be created dynamically as they are read.  If you want this, create
// them in a scope that is not shared by any [[BaseEnvironment]] object,
// and do the serialization through pointers yourself.
void BaseEnvironment::save()
{
  initialized=Environment::load;
  std::ofstream os(environment_file_fin,std::ios::binary);
  if (!os.is_open())
    throw Open_file_error(environment_file_fin);
  oarchive_t oa(os);

  scope_t& scope=scopes[scope_name];
  scope_t::size_type nenv=scope.size();
  oa << nenv;
  for (Environment* e : scope)
    e->vserial(oa);
}

void BaseEnvironment::load()
{
  std::ifstream is(environment_file_ini,std::ios::binary);
  if (!is.is_open())
    throw Open_file_error(environment_file_ini);
  iarchive_t iar(is);

  scope_t& scope=scopes[scope_name];
  scope_t::size_type nenv;
  iar >> nenv;
  if (nenv!=scope.size())
    throw Environment_unreadable("Number of objects has changed (now "+
				 std::to_string(scope.size()) + ", in file "+
   				 std::to_string(nenv)+")",HERE);
  for (Environment* e : scope)
    e->vserial(iar);
}

/******************************************************************************
 *
 * SimEnvironment
 *
 */

SimEnvironment_par::SimEnvironment_par(const char* scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("title",po::value<std::string>()->default_value("[untitled]"),"simulation title")
    ("max_steps",po::value<long>()->default_value(0),
     "maximum number of allowed steps in stage (0=unlimited)")
    ("log_interval",po::value<int>()->default_value(0),
     "interval (in steps) to write logs (0=no logging)")
    ;
}
  /// \name Public data initialized from parameters (SimEnvironment_par)
  /// @{
  std::string title;  ///< A title for the simulation
  long        max_steps;  ///< Maximum allowed steps in the stage (used by run()), 0=unlimited
  int         log_interval; ///<Interval for logging (in simulation steps)

  /// @} \name Public data computed/updated by Simulation::run() or step()
  /// @{
  long   steps_completed; ///< Steps completed since the first run
  long   steps_in_run;    ///< Steps completed in this run
  long   steps_in_stage;  ///< Steps completed in the present stage of this run
  bool   run_completed;   ///< Whether the run has been completed (set by Simulation::step())
  double time_completed;  ///< System time advanced since the first run, when applicable
  double time_in_run;     ///< System time advanced in this run, when applicable
  double time_in_stage;   ///< System time advanced in this stage, when applicable
  double precision;       ///< Precision reached, when applicable
  /// @}

void SimEnvironment::warm_init_local()
{
  BaseEnvironment::warm_init_local();
  title=par.value("title").as<std::string>();
  max_steps=par.value("max_steps").as<long>();
  log_interval=par.value("log_interval").as<int>();

  run_completed=false;
  steps_in_run=steps_in_stage=0;
  time_in_run=time_in_stage=0.;
}

void SimEnvironment::init_local()
{
  BaseEnvironment::init_local();
  title=par.value("title").as<std::string>();
  max_steps=par.value("max_steps").as<long>();
  log_interval=par.value("log_interval").as<int>();

  run_completed=false;
  steps_completed=steps_in_run=steps_in_stage=0;
  time_completed=time_in_run=time_in_stage=0.;
  precision=DBL_MAX;
}

}
