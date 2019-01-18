/*
 * stochastic.hh --  base environment for stochastic simulations
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

#ifndef STOCHASTIC_HH
#define STOCHASTIC_HH


// @ \chapter{Stochastic simulations}

// We provide parameters and environment useful as base for stochastic
// simulations.  Basically here we handle the random number generator.

#include <boost/serialization/split_member.hpp>
#include "parameters.hh"
#include "environment.hh"
#include "random.hh"

namespace glsim {

// @ \section{Stochastic parameters}

class StochasticParameters : public Parameters {
public:
  StochasticParameters(const char* scope=Parameters::default_scope);
} ;

inline StochasticParameters::StochasticParameters(const char* scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("stochastic.random_number_generator",po::value<std::string>()->required(),"generator to use")
    ("stochastic.seed",po::value<unsigned long int>()->default_value(0),
     "random generator seed")
    ;
}

// @ \section{Stochastic environment}

class StochasticEnvironment : public Environment {
public:
  StochasticEnvironment(const char* scope=Parameters::default_scope);
  virtual ~StochasticEnvironment();

  glsim::Random_number_generator *RNG;
  unsigned long seed;

protected:
  void init_local();
  void warm_init_local() {Environment::warm_init_local();}

private:
  StochasticParameters par;

  // Serialization
  friend class ::boost::serialization::access;
  template <typename Archive>
  void save(Archive &ar,const unsigned int version) const;
  template <typename Archive>
  void load(Archive &ar,const unsigned int version);
  BOOST_SERIALIZATION_SPLIT_MEMBER();

  virtual void vserial(iarchive_t &ar) {ar >> *this;}
  virtual void vserial(oarchive_t &ar) {ar << *this;}
} ;

// @ \paragraph{Construction and destruction}

inline StochasticEnvironment::StochasticEnvironment(const char* scope) :
  Environment(scope),
  RNG(0),
  seed(0),
  par(scope)
{}

inline StochasticEnvironment::~StochasticEnvironment()
{
  delete RNG;
}

template <typename Archive>
inline void StochasticEnvironment::save(Archive &ar,const unsigned int version) const
{
  ar << ::boost::serialization::base_object<Environment>(*this);
  ar << seed << RNG;
}

template <typename Archive>
inline void StochasticEnvironment::load(Archive &ar,const unsigned int version)
{
  ar >> ::boost::serialization::base_object<Environment>(*this);
  delete RNG;
  ar >> seed >> RNG;
}

} /* namespace */

BOOST_CLASS_VERSION(glsim::StochasticEnvironment,0);

#endif /* STOCHASTIC_HH */
