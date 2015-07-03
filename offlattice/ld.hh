/*
 * ld.hh -- Langevin dynamics simulations
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

#ifndef LD_HH
#define LD_HH

#include <unordered_map>

#include "mdenvironment.hh"
#include "stochastic.hh"
#include "simulation.hh"
#include "interactions.hh"

namespace glsim {

/*****************************************************************************/

class LDParameters : public Parameters {
public:
  LDParameters(const char* scope);
} ;

class LDEnvironment : public MDEnvironment {
public:
  LDEnvironment(const char* scope=Parameters::default_scope);

  ///@{ \name Set from parameters
  double  eta;
  double  temperature;
  
  ///@}@{ Computed by the simulation (guaranteed accurate only upon request)

  StochasticEnvironment SE;

protected:
  void init_local(),warm_init_local();
  
private:
  LDParameters par;

  void vserial(oarchive_t &ar) {ar << *this;}
  void vserial(iarchive_t &ar) {ar >> *this;}
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  friend class boost::serialization::access;

public:
  static const int class_version=1;
} ;

/** \class LDSimulation
    \ingroup OfflatticeSIM
    \brief Langevin dynamics simulation 

Sometimes called simply stochastic dynamics, this is a stochastic
dynamics simulation according to the Langevin equation with an inertia
term.

*/
class LDSimulation : public Simulation {
public:
  LDSimulation(LDEnvironment& e,OLconfiguration &c,Interactions *i);
  ~LDSimulation();
  const char  *name() const {return "Langevin Dynamics simulation";}
  void        step();
  void        log();
  void        log_start_sim();

protected:
  void update_observables();

  LDEnvironment&   env;
  OLconfiguration& conf;
  Interactions     *inter;

private:
  double Dt;
  std::unordered_map<short,double> c0,c1dt,c1mc2dt,c2dtsq,c2dt;
  std::unordered_map<short,BivariateGaussian_distribution*> noise;
} ;

} /* namespace */

BOOST_CLASS_VERSION(glsim::LDEnvironment,glsim::LDEnvironment::class_version);

#endif /* LD_HH */
