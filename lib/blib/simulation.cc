/*
 * simulation.cc -- definitions for simulation classes
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

#include <iomanip>
#include <sys/stat.h>
#include <unistd.h>
#include <netdb.h>
#include <arpa/inet.h>

#include "simulation.hh"

namespace glsim {

/******************************************************************************
 * 
 * Simulation
 * 
 */

// Construction and static data

volatile sig_atomic_t Simulation::termination_requested,
                      Simulation::signal_received;
struct sigaction Simulation::postpone_term;

/** 
We install a new handler for SIGTERM and SIGINT so that the sim can
terminate gracefully at the end of the current step.  The handler is
installed with the SA_RESTART flag so that library functions are
resumed after the handler ends.  Signals are blocked until handler is
installed.
*/
void Simulation::set_up_signals()
{
  sigemptyset(&postpone_term.sa_mask);
  sigaddset(&postpone_term.sa_mask,SIGTERM);
  sigaddset(&postpone_term.sa_mask,SIGINT);
  sigprocmask(SIG_BLOCK,&postpone_term.sa_mask,0);

  postpone_term.sa_handler=Simulation::sigterm_handler;
  postpone_term.sa_flags=SA_RESTART;
  sigaction(SIGTERM,&postpone_term,0);
  sigaction(SIGINT,&postpone_term,0);
  sigprocmask(SIG_UNBLOCK,&postpone_term.sa_mask,0);
}

void Simulation::sigterm_handler(int signal)
{
  sigprocmask(SIG_BLOCK,&postpone_term.sa_mask,0);
  termination_requested=GLSIM_TERM_ON_SIGNAL;
  signal_received=signal;
}

/**
run() executes the abstract simulation algorithm.

run() takes care of:

 - counts the number of steps completed in total, in the presen run
   and in the present stage of the run

 - calls step() repeatedly until either
   - a termination signal (SIGTERM or SIGINT) is received
   - the maximum number of steps (if specified) is reached
   - the maximum allowed wall time is reached (TO BE IMPLEMENTED)

 - calls the logging functions as appropriate

 - calls env.step() to advance the environment (also triggering
   observation of desired quantities)

 */
long Simulation::run()
{
  set_up_signals();
  log_start_sim();
  env.steps_in_stage=0;
  while (!env.run_completed && termination_requested==0) {
    env.steps_completed++;
    env.steps_in_run++;
    env.steps_in_stage++;
    step();
    env.step();
    if (env.steps_completed % env.log_interval==0) log();
    if (env.max_steps>0 && env.steps_in_stage>env.max_steps)
      termination_requested=GLSIM_TERM_ON_MAX_STEPS;
  }

  switch (termination_requested) {
  case 0:
    break;
  case GLSIM_TERM_ON_SIGNAL:
    logs(warn) << "\nWARNING: Terminating on signal " << signal_received
	       << "\n\n";
    break;
  case GLSIM_TERM_ON_MAX_STEPS:
    logs(warn) << "\nWARNING: Reached max_steps.\n\n";
      break;
  }
  log_stop_sim();
  return env.steps_in_stage;
}

/**
Prints and records starting time
*/
void Simulation::log_start_sim()
{
  char hn[1000];
  gethostname(hn,1000);
  struct hostent *he=gethostbyname(hn);
  logs(info) << "\nStarting simulation: " << name() << "\n\n";
  logs(info) << "Running on " << he->h_name << "[";
  for (struct in_addr **ap=(struct in_addr **)he->h_addr_list; *ap!=0; ++ap)
    logs(info) << " " << inet_ntoa(**ap);
  logs(info) << " ]\n";
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  boost::gregorian::date  d=now.date();
  logs(info) << "\n***** SIMULATION START ***** "
	     << d.day_of_week() << ' '
	     << d.day() << '-' << d.month() << '-' << d.year()
	     << ' ' << now.time_of_day()
	     << "\n\n";
  simulation_timer.start();
}

/**
   Computes elapsed time and writes to log stream
*/
void Simulation::log_stop_sim()
{
  using namespace boost::posix_time;
  using namespace boost::gregorian;
  using namespace boost::timer;

  simulation_timer.stop();
  ptime now = second_clock::local_time();
  date  d=now.date();
  logs(info) << "\n***** SIMULATION STOP ****** "
	     << d.day_of_week() << ' '
	     << d.day() << '-' << d.month() << '-' << d.year()
	     << ' ' << now.time_of_day()
	     << '\n';

  logs(info) << "\nCompleted " << env.steps_in_run << " steps in this run (" << env.steps_in_stage << 
      " in this stage)\n";

  cpu_times     simtimes=simulation_timer.elapsed();
  time_duration total=millisec(simtimes.wall/1000000);
  time_duration user=millisec(simtimes.user/1000000);
  time_duration system=millisec(simtimes.system/1000000);

  logs(info) << "\n      Elapsed time : "
	     << boost_duration_to_string(total) << '\n';
  logs(info) << "      CPU time     : "
	     << boost_duration_to_string(user+system)
	     << "  (" << boost_duration_to_string(user) << " user)\n";
  if (total.total_milliseconds()>0)
    logs(info) << "      Steps/second : "
	       << 1000.0*env.steps_in_stage / total.total_milliseconds() << "\n\n";
}

std::string 
Simulation::boost_duration_to_string(boost::posix_time::time_duration t)
{
  long days=t.hours()/24;
  long hours=t.hours() % 24;
  int  centisecs=100.0*t.fractional_seconds() /
    boost::posix_time::time_duration::ticks_per_second();
  std::stringstream ss;
  ss << days << ' ' << std::setfill('0') << std::setw(2) << hours << ':';
  ss << std::setw(2) << t.minutes() << ':' << std::setw(2) << t.seconds()
     << '.' << std::setw(2) << centisecs;
  return ss.str();
}

/******************************************************************************
 *
 * prepare
 *
 */


enum existing_file_status {none, partial, exist};

static existing_file_status check_existing_files(SimEnvironment &env,
						 bool ignore_partial);

void prepare(SimulationCL& CL,SimEnvironment &env,Configuration &conf)
{
  // Check for existing files and partial run
  env.init_base();
  existing_file_status old_files=
    check_existing_files(env,CL.value("ignore-partial-run").as<bool>());

  // Partial exist; load and resume
  if (old_files==partial) {
    logs(info) << "\nResuming previous run\n";
    conf.load(env.configuration_file_fin);  // env is already loaded
    return;
  }

  // Warn overwrite
  if (old_files==exist && CL.value("force-overwrite").as<bool>()==false)
    throw glsim::Runtime_error("Existing files would be overwritten, use -f to force");

  // Almost ready to go
  if (CL.value("initial_infix").as<std::string>()=="+++") {
    // This run starts from scratch; init everything
    env.init();
    CL.count("configuration-init")>0 ? 
      conf.init(CL.value("configuration-init").as<std::string>()) :
      conf.init();
  } else {
    // This is a continuation run; read and warm_init env, read conf
    env.load();
    env.warm_init();
    conf.load(env.configuration_file_ini);
  }
}

static existing_file_status check_existing_files(SimEnvironment &env,
						 bool ignore_partial)
{
  struct stat s;
  bool env_exist,conf_exist;

  logs(info) << "Checking for existing files...";
  // if (stat(envfile,&s)!=0) {
  if (stat(env.environment_file_fin.c_str(),&s)!=0) {
    logs(info) << " environment not found...";
    env_exist=false;
  } else {
    logs(info) << " environment found...";
    env_exist=true;
  }
  //  if (stat(conffile,&s)!=0) {
  if (stat(env.configuration_file_fin.c_str(),&s)!=0) {
    logs(info) << " configuration not found.\n";
    conf_exist=false;
  } else {
    logs(info) << " configuration found.\n";
    conf_exist=true;
  }

  if (!conf_exist && !env_exist) return none;   // Both missing
  if (!conf_exist || !env_exist) return exist;  // One of them missing, cannot continue

  if (ignore_partial) {
    logs(info) << "\n--ignore-partial-run given, not checking for partially completed run\n";
    return exist;
  }

  std::string fname=env.environment_file_ini;
  env.environment_file_ini=env.environment_file_fin;
  env.load();
  env.environment_file_ini=fname;
  if (env.run_completed) return exist;
  else {
    logs(info) << "Incomplete simulation found.\n";
    return partial;
  }
}


} /* namespace */
