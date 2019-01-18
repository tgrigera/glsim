/*
 * mdobservable.hh --  Recording basic MD quantities
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

#include <cstdio>

#include "olconfiguration.hh"
#include "glsim/observable.hh"
#include "md.hh"

#ifndef MDOBSERVABLE_HH
#define MDOBSERVABLE_HH

namespace glsim {

class MDObservable_parameters : public Parameters {
public:
  MDObservable_parameters(const char* scope);
} ;

class MDObservable : public SBObservable {
public:
  MDObservable(MDEnvironment&,OLconfiguration&);

  void interval_and_file();
  void write_header();
  void observe();

private:
  MDEnvironment           &env;
  OLconfiguration         &conf;
  MDObservable_parameters par;
} ;

inline MDObservable::MDObservable(MDEnvironment& e,OLconfiguration &c) :
  SBObservable(e),
  env(e),
  conf(c),
  par(e.scope())
{}

} /* namespace */

#endif /* MDOBSERVABLE_HH */
