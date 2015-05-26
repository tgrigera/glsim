/*
 * environment.hh -- declarations for environment classes
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

#ifndef ENVIRONMENT_HH
#define ENVIRONMENT_HH

#include <string>
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

The constructor must initialize the object in a suitable default, but
_must not_ attempt to read any value from Parameters: parameter
parsing will be done only after all environment objects are created.
Thus initialization from parameters must be provided with separate
init() methods, discussed further below.  Initialization comes in two
flavors: full initialization (init()), which rebuilds the environment
state from scratch, reading from the parameters object (and thus
file), and warm initialization, which means that the environment is
already initialized (typically, has been read from an earlier
simulation) but needs partial initialization to be prepared to start a
new simulation.  What exactly this means will depend on a simulation,
but for instance, in a Monte Carlo simulation warm initialization
would _not_ initialize the random number generator, while it _would_
initialize the requested number of steps.

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

*/
class Environment {
public:
  /// Create with default values and register in the scope
  Environment(const char* scope=Parameters::default_scope);
  virtual ~Environment();
  /// The scope we belong to
  const char * const scope() const {return scope_name.c_str();}

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
    throw Environment_wrong_version("CTSimEnvironment",version,class_version);
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
    throw Environment_wrong_version("CTSimEnvironment",version,class_version);
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

  std::string title;
  long        steps_completed,steps_in_run,steps_in_stage;
  long        max_steps;
  int         log_interval;
  bool        run_completed;

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
  static const int class_version=1;
} ;

inline SimEnvironment::SimEnvironment(const char* scope) :
  BaseEnvironment(scope),
  title("[untitled]"),
  steps_completed(0), steps_in_run(0),steps_in_stage(0),
  max_steps(0),
  log_interval(0),
  run_completed(false),
  par(scope)
{}

template <typename Archive>
inline void SimEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw Environment_wrong_version("CTSimEnvironment",version,class_version);
  ar & boost::serialization::base_object<BaseEnvironment>(*this);
  ar & title;
  ar & steps_completed & steps_in_run & steps_in_stage;
  ar & max_steps & log_interval;
  ar & run_completed;
}

/******************************************************************************
 *
 * CTSimEnviornment
 *
 */

/** \class CTSimEnvironment
    \ingroup Environment
    \brief Environment for continuous time simulations
*/
class CTSimEnvironment : public SimEnvironment {
public:
  CTSimEnvironment(const char* scope=Parameters::default_scope);

  double       time_completed,time_in_run,time_in_stage;

protected:
  void init_local();
  void warm_init_local();

private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  virtual void vserial(iarchive_t &ar) {ar >> *this;}
  virtual void vserial(oarchive_t &ar) {ar << *this;}

public:
  static const int class_version=1;
} ;

template <typename Archive>
inline void CTSimEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw Environment_wrong_version("CTSimEnvironment",version,class_version);
  ar & boost::serialization::base_object<SimEnvironment>(*this);
  ar & time_completed & time_in_run & time_in_stage;
}

} /* namespace */

BOOST_CLASS_VERSION(glsim::Environment,glsim::Environment::class_version);
BOOST_CLASS_VERSION(glsim::BaseEnvironment,glsim::BaseEnvironment::class_version);
BOOST_CLASS_VERSION(glsim::SimEnvironment,glsim::SimEnvironment::class_version);
BOOST_CLASS_VERSION(glsim::CTSimEnvironment,glsim::CTSimEnvironment::class_version);


#endif /* ENVIRONMENT_HH */
