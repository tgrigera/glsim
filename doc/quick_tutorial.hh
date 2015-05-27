/*
 * quick_tutorial.hh -- A quick tutorial on how to set up a simulation
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

/** \page Quick Quick tutorial

# A working simulation

To produce a working simulation based on `glsim` one must derive and
instantiate objects for the main abstractions of the library:

 - Parameters
 - Environment
 - Configuration
 - Simulation

Let's see how these objects must be defined.

## Configuration

Since this is a fake MD simulation, our configuration will be
OLconfiguration (which see).

~~~~~~~~~~{.cc}
#include "olconfiguration.hh"
~~~~~~~~~~

## Parameters and environment

The parameters classes represent the information that is read from a
.ini file.  They do not hold data, just declare the parameters
(i.e. the key=value pairs) that are accepted or expected.

Parameters needs to inherit from glsim::Parameters, and define the
parameters in its constructor.  Parsing will be handled
automatically if this is put in the main scope (the one that gets
sent to glsim::prepare()

~~~~~~~~~~{.cc}

#include "simulation.hh"

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
~~~~~~~~~~

The enviroment declares (typically public) data members wiht the
useful data it will hold.  The constructor initializes this to some
defaults.  A private par object of type MDPar is declared which is
initialized with the desired scope by the constructor.  Then the
following methods are needed:

- init_local() and warm_init_local():  Initalize from the par object (these methods will be called
                           when par is already properly initialized)

- the serialization functions: the two overloaded vserial functions (copy them verbatim), plus
                               serialize().  serialize MUST call base_serializa and then
                               serialize the class (see example below)

- the class version information (for serialization): the version
                                public static variable and the Boost macro

~~~~~~~~~~{.cc}
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

~~~~~~~~~~

## Simulation and observables

The simulation needs step() plus eventually overriding log

~~~~~~~~~~{.cc}

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
~~~~~~~~~~

## Main


 main() calls wmain through the convenience wrapper, then it must.

 - define objects for enviroment (all scopes) and configurations
 - call prepare (or interpret the command line)
 - initialize the scopes that were not passed to prepare()


~~~~~~~~~~{.cc}
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
  return glsim::UtilityEC(argc,argv,wmain);
}

~~~~~~~~~~

*/
