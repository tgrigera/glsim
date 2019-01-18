/*
 * environment.hh -- declarations for environment classes
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


#ifndef ENVIRONMENT_HH
#define ENVIRONMENT_HH

#include <string>
#include <list>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/version.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "log.hh"
#include "parameters.hh"

namespace glsim {
  
/** \class Environment
    \ingroup Environment

All environment types must be derived from this class.  This
implements the consolidation by scope name, and provides two functions
to build filenames according to the `glsim` convention.  Environment
is an abstract class because step() is a pure virtual.

The environment constructor must accept a scope name (a string), but
must provide a default name to allow default construction.  The scope
name will be also passed to the accompanying Parameters object.

The user of the library will derive at least one class from
BaseEnvironment below: that class has the full interface for loading
and saving complete scopes, and will be typically manged from the
simulation class or from the main program.  In principle, only one
object of this class per scope will be created.  Additionally and
optionally, one or more classes can be derived directly from
Environment: these will be the ``floating'' environments, which will
belong to the scope of one BaseEnvironment descendant, and saving,
loading and initialization will be managed from there.  For this
reason the initialization methods are kept protected at this level.

## Initialization

The constructor must initialize the object in a suitable default, but
_must not_ attempt to read any value from Parameters: parameter
parsing will be done only after all environment objects are created.

There are three ways initialization can happen in the simulation, and
the derived classes must be prepared to handle them correctly:

 1. _cold_ or _full_ initialization: This happens at the start of the
 first run, when there is no previously saved environment, and the
 environment must be completely initialized from the parameter files
 (.ini) reading through a Parameters object.  This is done by
 init_local() (which will be called by base BaseEnvironment::init).

 2. _warm_ or _partial_ initialization: This happen at the start of a
 continuation run, when the environment is already initialized
 (typically, has been read from an earlier simulation) but needs
 partial initialization to be prepared to start a new run that will be
 the logical dyamic continuation of the first.  What exactly this
 means will depend on a simulation, but for instance, in a Monte Carlo
 simulation warm initialization would _not_ initialize the random
 number generator, while it _would_ initialize the requested number of
 steps.  This must be done by warm_init_local() (called from
 BaseEnvironment::warm_init).

 3. initialization from a saved environment only: This happens when a
 partially completed run is found and a new simulation stage starts
 that will try to complete the previously started run.  No
 initialization method is involved, just BaseEnvironment::load, which
 loads each evironment through the serialize method.  Thus serialize
 must save and load all variables, even those that are read from
 Parameters and never changed.  This is necessary in order to avoid
 useless calls of init_local or warm_init_local, which may be
 expensive or have undesired side-effects.

 The protected variable initialized records the initialization method
 used and can be read through the public method initialization_kind().

*/
class Environment {
public:
  /// Create with default values and register in the scope
  Environment(const char* scope=Parameters::default_scope);
  virtual ~Environment();
  /// The scope we belong to
  const char * const scope() const {return scope_name.c_str();}

  enum init_method {deflt,warm,cold,load};
  init_method initialization_kind() const {return initialized;}
  

protected:
  typedef boost::archive::binary_oarchive oarchive_t;
  typedef boost::archive::binary_iarchive iarchive_t;
  // typedef boost::archive::text_oarchive oarchive_t;
  // typedef boost::archive::text_iarchive iarchive_t;

  /// Evolve the environment (noop at this level)
  virtual void step_local() {}
  /// Initialize from parameters file (just this object)
  virtual void init_local();
  /// Warm initialize from parameters file (just this object)
  virtual void warm_init_local();

  /// Build initial filename from arguments
  void initial_filename(std::string& name,const std::string& prefix);
  /// Build final filename from arguments
  void final_filename(std::string& name,const std::string& prefix);

  init_method initialized;  /// Initialization method used
  std::string scope_name;
  std::string ini_infix,fin_infix,extension;

private:
  Parameters par;

  typedef std::list<Environment*> scope_t;
  typedef std::map<std::string,scope_t> scopemap_t;
  static  scopemap_t scopes;
  
  friend class BaseEnvironment;
  virtual void vserial(iarchive_t &ar);
  virtual void vserial(oarchive_t &ar);
  friend class boost::serialization::access;

  void    init_common();
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);

public:
  static const int class_version=1;
} ;
  
/** \ingroup Exception
    \brief Thrown by Environment when unable to read from file
*/
class Environment_unreadable : public glsim::Runtime_error {
public:
  explicit Environment_unreadable(const std::string& msg,
				  const Source_context &c=Source_context()) :
      Runtime_error("ERROR: cannot read environmet: "+msg,c) {}
  ~Environment_unreadable() throw() {}
} ;

/** \ingroup Exception
    \brief Thrown by Environment when current version is incompatible with data stored in file
*/
class Environment_wrong_version : public glsim::Environment_unreadable {
public:
  explicit Environment_wrong_version(const std::string& classname,
				  const int version_in_file,
				  const int current_version,
				  const Source_context &c=Source_context()) :
      Environment_unreadable("incompatible version of class "+classname+
			     " (in file "+::std::to_string(version_in_file)+
			     " current "+::std::to_string(current_version)+")",c) {}
  ~Environment_wrong_version() throw() {}
} ;

/**

Serialization is done with `Boost::serialization`, thus a serialize()
function is needed.  However, to implement loading and saving at the
scope level, two other functions must be provided, although they are
one-liners.  The reason is that when saving a scope the serialize
functions will be called on references, which are not handled by Boost
(virtual functions are needed for serializing objects from a base
class reference, but the mechanism used by Boost employs template
functions, which cannot be virtual).  The functions vserial() below
will simply call the correct serialize function, but **must** be
overridden (though written identically), otherwise only the base part
will be saved (see tutorial).
*/
template <typename Archive>
inline void Environment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw Environment_wrong_version("Environment",version,class_version);
  ar & initialized;
  ar & scope_name;
  ar & ini_infix;
  ar & fin_infix;
}

inline void Environment::vserial(oarchive_t &ar) {ar << *this;}
inline void Environment::vserial(iarchive_t &ar) {ar >> *this;}

/******************************************************************************
 *
 * BaseEnvironment
 *
 */

/** \class BaseEnv_par
    \ingroup Environment
    \brief Parameters for BaseEnvironment
*/
class BaseEnv_par : public Parameters {
public:
  BaseEnv_par(const char *scope=Parameters::default_scope);
} ;

/** \class BaseEnvironment
    \ingroup Environment

This class is the base for simulation environments, i.e. environments
that can do loading, saving and initialization at the scope level.  We
thus store here the filenames for reading and writing, and implement
the scope-wide methods.  The local step method is empty, the
scope-wide step() calls step_local() for all members of the scope.
*/
class BaseEnvironment : public Environment {
public:
  BaseEnvironment(const char* scope=Parameters::default_scope);

  void parse(const char*filename);
  void init_base();
  /// Initialize whole scope from parameters file
  void init();
  /// Warm initialize whole scope from parameters file
  void warm_init();
  /// Scope-wide advance one step
  void step();
  /// Load whole scope from file
  void load();
  /// Write whole scope to file
  void save();

  /// Initial and final configuration files
  std::string configuration_file_ini,configuration_file_fin;
  /// Initial and final environmnent files
  std::string environment_file_ini,environment_file_fin;

protected:
  void init_local();
  void warm_init_local();

private:
  std::string configuration_file_prefix,environment_file_prefix;
  BaseEnv_par par;

  void init_local_common();

  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  virtual void vserial(iarchive_t &ar) {ar >> *this;}
  virtual void vserial(oarchive_t &ar) {ar << *this;}

public:
  static const int class_version=1;
} ;


/**
The constructor loads everything with suitable defaults.  For
initialization, we need first to write the corresponding parameters
constructor, which declares the parameters to Boost::program_options.
*/
inline BaseEnvironment::BaseEnvironment(const char* scope) :
  Environment(scope),
  configuration_file_ini("[AUTO]"),
  configuration_file_fin("[AUTO]"),
  environment_file_ini("[AUTO]"),
  environment_file_fin("[AUTO]"),
  configuration_file_prefix("conf"),
  environment_file_prefix("env"),
  par(scope)
{
}

/// Serialization
template <typename Archive>
inline void BaseEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw Environment_wrong_version("BaseEnvironment",version,class_version);
  ar & boost::serialization::base_object<Environment>(*this);
  ar & configuration_file_ini & configuration_file_fin;
  ar & environment_file_ini & environment_file_fin;
  ar & configuration_file_prefix & environment_file_prefix;
}

/******************************************************************************
 *
 * SimEnvironment
 *
 */
/** \class SimEnvironment_par
    \ingroup Environment
    \brief Parameters for SimEnvironment
*/
class SimEnvironment_par : public Parameters {
public:
  SimEnvironment_par(const char* scope=Parameters::default_scope);
} ;

/** \class SimEnvironment
    \ingroup Environment

This is the environment that needed by Simulation.  We store here
variables to keep track of the number of steps performed (total, in
the run, in the stage), and the desired log interval (log_interval).
There is also the boolean run_completed to communicate with
Simulation::run(), which is however not updated here, but in one of
the specialized environments or in Simulation::step().

The local step function is a no-op., since count of the number of
steps is kept in Simulation::run().  This can be more efficient,
because BaseEnvironment::step() can be told not call no-ops at all
(TODO). This function is typically overridden anyway.
*/
class SimEnvironment : public BaseEnvironment {
public:
  SimEnvironment(const char* scope=Parameters::default_scope);

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

protected:
  void init_local();
  void warm_init_local();

private:
  SimEnvironment_par par;

  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  virtual void vserial(iarchive_t &ar) {ar >> *this;}
  virtual void vserial(oarchive_t &ar) {ar << *this;}

public:
  static const int class_version=2;
} ;

inline SimEnvironment::SimEnvironment(const char* scope) :
  BaseEnvironment(scope),
  title("[untitled]"),
  max_steps(0),
  log_interval(0),
  steps_completed(0), steps_in_run(0), steps_in_stage(0),
  run_completed(false),
  time_completed(0.), time_in_run(0.), time_in_stage(0.),
  precision(0.),
  par(scope)
{}

template <typename Archive>
inline void SimEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw Environment_wrong_version("SimEnvironment",version,class_version);
  ar & boost::serialization::base_object<BaseEnvironment>(*this);
  ar & title;
  ar & max_steps & log_interval;
  ar & steps_completed & steps_in_run & steps_in_stage;
  ar & run_completed;
  ar & time_completed & time_in_run & time_in_stage;
  ar & precision;
}

} /* namespace */

BOOST_CLASS_VERSION(glsim::Environment,glsim::Environment::class_version);
BOOST_CLASS_VERSION(glsim::BaseEnvironment,glsim::BaseEnvironment::class_version);
BOOST_CLASS_VERSION(glsim::SimEnvironment,glsim::SimEnvironment::class_version);

#endif /* ENVIRONMENT_HH */
