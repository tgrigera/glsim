/*
 * observable.cc --  base class for observable quantities (definitions)
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

#include "observable.hh"

namespace glsim {

void SBObservable::step_local()
{
  if (obs_interval==0) return;
  if (of==0) {
    of=fopen(obs_fname.c_str(),"a");
    if (!of) throw Open_file_error(obs_fname);
  }
  if (senv.steps_completed % obs_interval==0) observe();
}

void KMCObservable::init_local()
{
  Observable<double>::init_local();
  obs_time=0.;
}

void KMCObservable::step_local()
{
  if (obs_interval==0) return;
  if (of==0) {
    of=fopen(obs_fname.c_str(),"a");
    if (!of) throw Open_file_error(obs_fname);
  }
  double times=(env.time_completed-obs_time)/obs_interval;
  if (times>1000) {
     fprintf(of,"WARNING Sthg not good!!!\ntoo many observations (%g), time %g\n",times,env.time_completed);
     return;
  }
  while (env.time_completed > obs_time) {
    observe();
    obs_time+=obs_interval;
  }
}


} /* namespace */
