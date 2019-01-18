/*
 * trajectory.cc --  observer for saving trajectories along a run
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

#include "trajectory.hh"

namespace glsim {

void Trajectory::init_local()
{
  Environment::init_local();
  first_observed=false;
  common_init();
}

void Trajectory::warm_init_local()
{
  Environment::warm_init_local();
  common_init();
}

void Trajectory::common_init()
{
  traj_fname="[AUTO]";
  obs_interval=par.value("trajectory.obs_interval").as<int>();
  if (obs_interval>0) {
    traj_file_prefix=par.value("trajectory.file_prefix").as<std::string>();
    final_filename(traj_fname,traj_file_prefix);
  }
}

/// This *must* be called before the first simulation step (but after
/// simulation construction).
void Trajectory::observe_first()
{
  if (obs_interval==0) return;
  if (!of) {
    iconf=new OLconfiguration();
    iconf->N=conf.N;
    file_opt.time_frame();
    memcpy(iconf->box_length,conf.box_length,3*sizeof(double));
    memcpy(iconf->box_angles,conf.box_angles,3*sizeof(double));
    iconf->id=conf.id;
    if (par.value("trajectory.record_type").as<bool>())
      iconf->type=conf.type;
    if (par.value("trajectory.record_flags").as<bool>())
      iconf->flags=conf.flags;
    if (par.value("trajectory.record_r").as<bool>()) {
      iconf->r=conf.r;
      file_opt.r_frame();
    }
    if (par.value("trajectory.record_v").as<bool>()) {
      iconf->v=conf.v;
      file_opt.v_frame();
    }
    if (par.value("trajectory.record_a").as<bool>()) {
      iconf->a=conf.a;
      file_opt.a_frame();
    }
    of=new OLconfig_file(iconf,file_opt);
    of->create(traj_fname.c_str());
    of->write_header();
  }
  if (first_observed) return;
  step_local();
  first_observed=true;
}

void Trajectory::step_local()
{
  if (obs_interval==0) return;
  if (senv.steps_completed % obs_interval==0) {
    iconf->time=conf.time;
    iconf->step=conf.step;
    memcpy(iconf->box_length,conf.box_length,3*sizeof(double));
    memcpy(iconf->box_angles,conf.box_angles,3*sizeof(double));
    of->append_record();
  }
}

} /* namespace */
