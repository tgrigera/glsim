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


namespace glsim {

#define GLSIM_TERM_ON_SIGNAL    1
#define GLSIM_TERM_ON_MAX_STEPS 2

/** \class Simulation
    \ingroup Simulation

This is the base for all (nonparallel at least) simulations.  It
handles logging and termination (by interruption through signals or on
completion of the simulation).
*/
class Simulation {
public:
  /// The constructor just stores references to environment and configuration
  Simulation(glsim::SimEnvironment&,glsim::Configuration&);
  virtual const char *name() const=0;
  /// This runs the simulation, calling step() as appropriate
  virtual long run();

  /// Step must be defined by the actual simulation class
  virtual void step()=0;
  /// Log start of simulation and record initial time
  virtual void log_start_sim();
  /// Log end of simulation and print elapsed and final time
  virtual void log_stop_sim();
  virtual void log() {}
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

/** \ingroup Simulation
    \brief Preparing configuration and environment

Configuration and environment must be ready to run before the
simulation is created.  The following function can be used to set them
up.  They rely on SimulationCL to find the appropriate files from the
command line and through partial initialization of the environment.
Before the prepare() call, all the environments belonging to the scope
of the environment passed to prepare() must have been constructed in
their default state, or the call will fail (perhaps silently and
miserably).

 \param CL  The command-line parameters (must be already parsed)
 \param env The eviroment to be initialized
 \param conf The cofiguration to be initialized

*/
void prepare(SimulationCL &CL,SimEnvironment &env,Configuration &conf);


} /* namespace */

BOOST_CLASS_VERSION(glsim::Simulation,0);

#endif /* SIMULATION_HH */
