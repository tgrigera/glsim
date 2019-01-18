/*
 * mcobservable.hh --  Recording basic quantities along a MC run
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

#ifndef MCOBSERVABLE_HH
#define MCOBSERVABLE_HH

#include <cstdio>

#include "olconfiguration.hh"
#include "glsim/observable.hh"
#include "mc.hh"

namespace glsim {

class MCObservable_parameters : public Parameters {
public:
  MCObservable_parameters(const char* scope);
} ;

class MCObservable : public SBObservable {
public:
  MCObservable(MCEnvironment&,OLconfiguration&);

  void interval_and_file();
  void write_header();
  void observe();

private:
  MCEnvironment           &env;
  OLconfiguration         &conf;
  MCObservable_parameters par;
} ;

inline MCObservable::MCObservable(MCEnvironment& e,OLconfiguration &c) :
  SBObservable(e),
  env(e),
  conf(c),
  par(e.scope())
{}

} /* namespace */

#endif /* MCOBSERVABLE_HH */
