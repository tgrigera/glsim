/*
 * mdenvironment.hh -- Environment for Molecular Dynamics
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

#ifndef MDENVIRONMENT_HH
#define MDENVIRONMENT_HH

#include "environment.hh"

namespace glsim {

class MDParameters : public Parameters {
public:
  MDParameters(const char *scope);
} ;

/** \class MDEnvironment
    \ingroup OfflatticeSIM
*/
class MDEnvironment : public SimEnvironment {
public:
  MDEnvironment(const char* scope=Parameters::default_scope);

  long   MDsteps;
  double time_step;

protected:
  void init_local(),warm_init_local();
  
private:
  MDParameters par;

  void vserial(oarchive_t &ar) {ar << *this;}
  void vserial(iarchive_t &ar) {ar >> *this;}
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  friend class boost::serialization::access;

public:
  static const int class_version=1;
} ;

} /* namespace */


BOOST_CLASS_VERSION(glsim::MDEnvironment,glsim::MDEnvironment::class_version);

#endif /* MDENVIRONMENT_HH */
