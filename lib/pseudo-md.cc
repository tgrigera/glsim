/*
 * pseudo-md.cc -- a fake MD simulation for tutorial purposes
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

#include "simulation.hh"

/// Configuration is OLConfiguration
#include "olconfiguration.hh"

/******************************************************************************
 * Parameters and Environment
 *
 */

// Parameters needs to inherit from glsim::Parameters, and define the
// parameters in its constructor.  Parsing will be handled
// automatically if this is put in the main scope (the one that gets
// sent to glsim::prepare()
class MDPar : public glsim::Parameters {
public:
  MDPar(const char *scope);
} ;

MDPar::MDPar(const char* s) :
  Parameters(s)
{
  parameter_file_options().add_options()
    ("MD.time_step",po::value<double>()->required(),"time step for MD integrator")
    ("MD.steps",po::value<long>()->required(),"number of steps to run")
    ;
}

// The enviroment declares (typically public) data members wiht the
// useful data it will hold.  The constructor initializes this to some
// defaults.  A private par object of type MDPar is declared which is
// initialized with the desired scope by the constructor.  Then the
// following methods are needed:
//
// - init_local() and warm_init_local():  Initalize from the par object (these methods will be called
//                            when par is already properly initialized)
//
// - the serialization functions: the two overloaded vserial functions (copy them verbatim), plus
//                                serialize().  serialize MUST call base_serializa and then
//                                serialize the class (see example below)
//
// - the class version information (for serialization): the version
//                                 public static variable and the Boost macro
//
class MDEnvironment : public glsim::SimEnvironment {
public:
  MDEnvironment(const char* scope=glsim::Parameters::default_scope);

  double time_step;
  long   MDsteps;

protected:
  void init_local(),warm_init_local();

private:
  MDPar par;

  void vserial(oarchive_t &ar) {ar << *this;}
  void vserial(iarchive_t &ar) {ar >> *this;}
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  friend class boost::serialization::access;

public:
  static const int class_version=2;
} ;

inline MDEnvironment::MDEnvironment(const char* scope) :
  SimEnvironment(scope),
  time_step(0),
  MDsteps(0),
  par(scope)
{}

// DO NOT forget to call the base's init!!
void MDEnvironment::init_local()
{
  SimEnvironment::init_local();
  time_step=par.value("MD.time_step").as<double>();
  MDsteps=par.value("MD.steps").as<long>();
}

void MDEnvironment::warm_init_local()
{
  SimEnvironment::warm_init_local();
  time_step=par.value("MD.time_step").as<double>();
  MDsteps=par.value("MD.steps").as<long>();
}

// DO NOT forget the base_object call!
template <typename Archive>
inline void MDEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw glsim::Environment_wrong_version("MDEnvironment",version,class_version);
  ar & boost::serialization::base_object<SimEnvironment>(*this);
  ar & time_step;
  ar & MDsteps;
}

BOOST_CLASS_VERSION(MDEnvironment,MDEnvironment::class_version);

/******************************************************************************
 * 
 * Simulation and observables
 *
 */

/*

  The simulation needs step() plus eventually overriding log

 */

class MDSim : public glsim::Simulation {
public:
  MDSim(MDEnvironment& e,glsim::OLconfiguration &c) :
    Simulation(e,c),
    env(e), conf(c)
  {}
  const char* name() const { return "Fake MD";}
  void step();

private:
  MDEnvironment& env;
  glsim::OLconfiguration& conf;
} ;


//
// step must set env.run_completed if we are finished (or have someone else do it
//
// note that env.steps_completed etc are actually one step ahead on
// entering step (i.e. they hold the number of steps that will be
// current at the end of the present call to step).

void MDSim::step()
{
  env.time_completed+=env.time_step;
  env.time_in_run+=env.time_step;
  env.run_completed = env.steps_in_run>=env.MDsteps;

  sleep(1);

  conf.step=env.steps_completed;
  conf.time=env.time_completed;
  glsim::logs(glsim::debug) << "Steps " << conf.step << " time " << env.time_completed << '\n';


}

/******************************************************************************
 * 
 * main
 *
 */

/*

  main() calls wmain through the convenience wrapper, then it must.

  - define objects for enviroment (all scopes) and configurations
  - call prepare (or interpret the command line)
  - initialize the scopes that were not passed to prepare()

 */
void wmain(int argc, char *argv[])
{
  MDEnvironment  env;
  glsim::OLconfiguration conf;
  glsim::SimulationCL CL("pseudo-MD","(C) 2015 Tomas S. Grigera",env.scope());
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  // Initialize other scopes if they exist
  MDSim sim(env,conf);
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
