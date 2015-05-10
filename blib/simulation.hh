/*
 * simulation.hh -- declarations for simulation classes
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

#ifndef SIMULATION_HH
#define SIMULATION_HH

#include <csignal>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_duration.hpp>
#include <boost/timer/timer.hpp>

#include "log.hh"
#include "environment.hh"
#include "configuration.hh"


// @ \chapter{Simulation}

// - sim::run() executes the abstract simulation algorithm

// Define step-based and target-based simulations; explain how both can
// be dealt with here.  env.run_completed must be set by child sim or
// specific environment.

// Before the simulation is created, Env and configuration must be ready
// to run.   Simulation wil \emph{not} initialize config or env.  To aid
// in this initialization we provide a [[prepare]] function below.


namespace glsim {

#define GLSIM_TERM_ON_SIGNAL    1
#define GLSIM_TERM_ON_MAX_STEPS 2

class Simulation {
public:
  Simulation(glsim::SimEnvironment&,glsim::Configuration&);
  virtual const char *name() const=0;
  virtual long run();

  virtual void step()=0;
  virtual void log_start_sim();
  virtual void log() {}
  virtual void log_stop_sim();
  virtual ~Simulation() {}

protected:
  void   set_up_signals();

  SimEnvironment &env;
  Configuration  &conf;

  static volatile sig_atomic_t termination_requested,signal_received;

private:
  static void sigterm_handler(int);
  std::string boost_duration_to_string(boost::posix_time::time_duration);
  static struct sigaction postpone_term;
  boost::timer::cpu_timer simulation_timer;
} ;

inline Simulation::Simulation(glsim::SimEnvironment& e,glsim::Configuration& c) :
  env(e),conf(c)
{}

// @ \section{Preparing configuration and environment}

// Configuration and environment must be ready to run before the
// simulation is created.  The following function can be used to set them
// up.  They rely on [[StandardCL]] to find the appropriate files from
// the command line and through partial initialization of the
// environment.  Before the [[prepare()]] call, all the environments
// belonging to the scope of the environment passed to [[prepare()]] must
// have been constructed in their default state, or the call will fail
// (perhaps silently and miserably).
void prepare(int argc,char *argv[],SimEnvironment &env,Configuration &conf);


} /* namespace */

BOOST_CLASS_VERSION(glsim::Simulation,0);

#endif /* SIMULATION_HH */
