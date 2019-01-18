/*
 * mdenvironment.hh -- Environment for Molecular Dynamics
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

#ifndef MDENVIRONMENT_HH
#define MDENVIRONMENT_HH

#include "environment.hh"
#include "interactions.hh"

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

  ///@{ \name Set from parameters
  long   MDsteps;
  double time_step;

  ///@}@{ Computed by the simulation (guaranteed accurate only upon request)
  double      Etot;    ///< Total energy
  double      Ekin;    ///< Kinetic energy
  double      Epot;    ///< Potential energy
  double      Ptot[3]; ///< Total momentum
  double      total_mass;
  double      total_number;

  Interactions *inter;  ///< The simulation can/should place here the pointer to interactions to allow observers to report interactions information
  ///@}

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
