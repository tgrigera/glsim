/*
 * observable.cc --  base class for observable quantities (definitions)
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
