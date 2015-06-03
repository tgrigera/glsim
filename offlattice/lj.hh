/*
 * lj.hh -- The Lennard-Jones potential
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

#ifndef LJ_JH
#define LJ_JH

#include <math.h>

#include "environment.hh"

namespace glsim {

/*****************************************************************************
 *
 * Repulsive Lennard-Jones
 *
 */

//
// Parameters
//

class RepulsiveLennardJones_parameters : public Parameters {
public:
  RepulsiveLennardJones_parameters(const char* scope);
} ;

inline RepulsiveLennardJones_parameters::
RepulsiveLennardJones_parameters(const char* scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("RLJ.sigma",po::value<double>(),"LJ length scale")
    ("RLJ.epsilon",po::value<double>(),"LJ energy scale")
    ("RLJ.mass",po::value<double>(),"particle mass")
    ;
}

//
// Environment
//

class RepulsiveLennardJones_environment : public Environment {
public:
  RepulsiveLennardJones_environment(const char* scope=glsim::Parameters::default_scope);

  double sigma,epsilon,mass;

protected:
  void init_local(),warm_init_local();

private:
  RepulsiveLennardJones_parameters par;
  
  void vserial(oarchive_t &ar) {ar << *this;}
  void vserial(iarchive_t &ar) {ar >> *this;}
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  friend class boost::serialization::access;

public:
  static const int class_version=1;
} ;

inline RepulsiveLennardJones_environment::
RepulsiveLennardJones_environment(const char* scope) :
  Environment(scope),
  sigma(1.),
  epsilon(1.),
  mass(1.),
  par(scope)
{}

inline void RepulsiveLennardJones_environment::init_local()
{
  Environment::init_local();
  sigma=par.value("RLJ.sigma").as<double>();
  epsilon=par.value("RLJ.epsilon").as<double>();
  mass=par.value("RLJ.mass").as<double>();
}

inline void RepulsiveLennardJones_environment::warm_init_local()
{
  Environment::warm_init_local();
}

template <typename Archive>
void RepulsiveLennardJones_environment::
serialize(Archive &ar,const unsigned int version)
{
  ar & boost::serialization::base_object<Environment>(*this);
  ar & sigma & epsilon & mass;
}


/** \class RepulsiveLennardJones
    \ingroup OfflatticeINT

    This is the Lennard-Jones potential without the attractive tail
    (it is cut-off at the minimum).
*/
class RepulsiveLennardJones {
public:
  RepulsiveLennardJones(const char *scope);
  void init(OLconfiguration &c);
  const char* name() const {return "Purely repulsive Lennard-Jones particles";}
  bool within_cutoff(double dsq,int ta,int tb)
  {return dsq<=rcsq;}
  double mass(short t) {return env.mass;}
  double cutoff() {return rc;}

  double external_field(double *d,short t) {return 0;}
  double external_field(double *d,short t,double *f) { f[0]=f[1]=f[2]=0.; return 0.;}
  double pair_potential(double dsq,short t1,short t2);
  double pair_potential_r(double dsq,short t1,short t2);
  double pair_potential(double ds1,short t1,short t2,double &vpsr);
  double pair_potential_r(double ds1,short t1,short t2,double &vpsr);

private:
  RepulsiveLennardJones_environment env;

  double sigmasq,foureps,derivpf,rc,rcsq;
} ;

inline RepulsiveLennardJones::
RepulsiveLennardJones(const char *scope) :
  env(scope)
{
  sigmasq=env.sigma*env.sigma;
  foureps=4.*env.epsilon;
  derivpf=48.*env.epsilon/sigmasq;
  rc=pow(2,1./6)*env.sigma;
  rcsq=rc*rc;
}

inline void RepulsiveLennardJones::
init(OLconfiguration &c)
{
  sigmasq=env.sigma*env.sigma;
  foureps=4.*env.epsilon;
  derivpf=48.*env.epsilon/sigmasq;
  rc=pow(2,1./6)*env.sigma;
  rcsq=rc*rc;

  logs(info) << "Lennard-Jones parameters from environment are:\n"
	     << "  epsilon = " << env.epsilon
	     << "  sigma   = " << env.sigma
	     << "  mass    = " << env.mass
	     << "\n [applying cut-off at sigma^(1/6)]\n";

  // check configuration here

  if (c.type==0) {
    c.type=new short[c.N];
    memset(c.type,0,c.N*sizeof(short));
  }
}

inline
double RepulsiveLennardJones::pair_potential(double rsq,short type1,short type2)
{
  return rsq<rcsq ? pair_potential_r(rsq,type1,type2) : 0;
}

inline double RepulsiveLennardJones::
pair_potential_r(double rsq,short type1,short type2)
{
  double s12,s6;
  s6=sigmasq/rsq;
  s6*=s6*s6;
  s12=s6*s6;
  return foureps*(s12-s6)+env.epsilon;
}

/**
   This must return v'(r)/r
*/
inline double RepulsiveLennardJones::
pair_potential(double rsq,short t1,short t2,double &vpsr)
{
  if (rsq<rcsq)
    return pair_potential_r(rsq,t1,t2,vpsr);
  else {
    vpsr=0;
    return 0;
  }
}

inline double RepulsiveLennardJones::
pair_potential_r(double rsq,short t1,short t2,double &vpsr)
{
  double srsq,s12,s6;
  srsq=sigmasq/rsq;
  s6=srsq*srsq*srsq;
  s12=s6*s6;
  vpsr=derivpf*srsq*(s6/2.-s12);
  return foureps*(s12-s6)+env.epsilon;
}

} /* namespace */

BOOST_CLASS_VERSION(glsim::RepulsiveLennardJones_environment,
		    glsim::RepulsiveLennardJones_environment::class_version);

#endif /* LJ_JH */
