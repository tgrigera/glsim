/*
 * trajectory.hh --  observer for saving trajectories along a run
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

#ifndef TRAJECTORY_HH
#define TRAJECTORY_HH

#include "environment.hh"

#include "olconfiguration.hh"

namespace glsim {

class Trajectory_parameters : public Parameters {
public:
  Trajectory_parameters(const char* scope=Parameters::default_scope) :
    Parameters(scope)
  {parameter_file_options().add_options()
      ("trajectory.obs_interval",po::value<int>()->default_value(0),
       "observation interval, 0=don't observe")
      ("trajectory.file_prefix",po::value<std::string>(),
       "file prefix for trajectories")
      ("trajectory.record_type",po::bool_switch()->default_value(false),
       "whether to write type to the trajectory file")
      ("trajectory.record_flags",po::bool_switch()->default_value(false),
       "whether to write flags to the trajectory file")
      ("trajectory.record_r",po::bool_switch()->default_value(false),
       "whether to write positions to the trajectory file")
      ("trajectory.record_v",po::bool_switch()->default_value(false),
       "whether to write velocities to the trajectory file")
      ("trajectory.record_a",po::bool_switch()->default_value(false),
       "whether to write accelerations to the trajectory file")
      ;
  }
} ;

/** \class Trajectory
    \ingroup Observable

WARNING: This class *needs* a call to observe_first() to complete
initialization.

N, time, step, box(*2) and id are always written (and made frame if
asked on construction).  type flags r,v,a must be requested from
parameters.

*/
class Trajectory : public Environment {
public:
  Trajectory(const SimEnvironment&,const OLconfiguration&,
	     OLconfig_file::options);
  virtual ~Trajectory();
  void observe_first();
  void step_local();

protected:
  void init_local();
  void warm_init_local();

  OLconfig_file *of;
  std::string   traj_file_prefix,traj_fname;
  int           obs_interval;  

private:
  bool                   first_observed;
  const OLconfiguration  &conf;
  OLconfiguration        *iconf;
  OLconfig_file::options file_opt;
  const SimEnvironment&  senv;
  Trajectory_parameters  par;

  void                  common_init();
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  virtual void vserial(iarchive_t &ar) {ar >> *this;}
  virtual void vserial(oarchive_t &ar) {ar << *this;}
} ;

inline Trajectory::Trajectory(const SimEnvironment &e,const OLconfiguration& c,
		       OLconfig_file::options opt) :
  Environment(e.scope()),
  of(0),
  traj_file_prefix("traj"),
  traj_fname("[AUTO]"),
  obs_interval(0),
  first_observed(false),
  conf(c),
  iconf(0),
  file_opt(opt),
  senv(e),
  par(e.scope())
{}

inline Trajectory::~Trajectory()
{
  delete of;
  if (iconf) {
    iconf->id=0;
    iconf->type=0;  // Because these buffers are conf's, not ours
    iconf->flags=0;
    iconf->r=0;
    iconf->v=0;
    iconf->a=0;
    delete iconf;
  }
}

template <typename Arch>
void Trajectory::serialize(Arch &ar, const unsigned version)
{
  ar & traj_fname & obs_interval;
  ar & first_observed;
}

} /* namespace */

#endif /* TRAJECTORY_HH */
