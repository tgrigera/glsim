/*
 * observable.hh --  base class for observable quantities (declarations)
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

#ifndef OBSERVABLE_HH
#define OBSERVABLE_HH

#include "environment.hh"

namespace glsim {

// step_local() missing at this point.

/** \class Observable
    \brief Base for classes representing observables
    \ingroup Observable
*/
template <typename ITYP>
class Observable : public Environment {
public:
  Observable(SimEnvironment &);
  virtual ~Observable();
  void observe_first();

protected:
  void init_local();
  void warm_init_local();

  virtual void observe()=0;
  virtual void write_header()=0;
  virtual void interval_and_file()=0;

  FILE*       of;
  std::string obs_file_prefix,obs_fname;
  ITYP        obs_interval;  

private:
  SimEnvironment   &senv;
  bool             first_observed;

  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  virtual void vserial(iarchive_t &ar) {ar >> *this;}
  virtual void vserial(oarchive_t &ar) {ar << *this;}
} ;

template <typename ITYP>
Observable<ITYP>::Observable(SimEnvironment &e) :
  Environment(e.scope()),
  of(0),
  obs_file_prefix("obs"),
  obs_fname("[AUTO]"),
  obs_interval(0),
  senv(e),
  first_observed(false)
{}

template <typename ITYP>
Observable<ITYP>::~Observable()
{
  if (of) fclose(of);
}

template <typename ITYP>
template <typename Arch>
void Observable<ITYP>::serialize(Arch &ar, const unsigned version)
{
  ar & obs_fname & obs_interval;
  ar & first_observed;
}

template <typename ITYP>
void Observable<ITYP>::init_local()
{
  Environment::init_local();
  interval_and_file();
  obs_fname="[AUTO]";
  final_filename(obs_fname,obs_file_prefix);
  of=fopen(obs_fname.c_str(),"w");
  if (!of) throw Open_file_error(obs_fname);
  first_observed=false;
  write_header();
}

template <typename ITYP>
void Observable<ITYP>::warm_init_local()
{
  Environment::warm_init_local();
  interval_and_file();
  obs_fname="[AUTO]";
  final_filename(obs_fname,obs_file_prefix);
  of=fopen(obs_fname.c_str(),"w");
  if (!of) throw Open_file_error(obs_fname);
  first_observed=true;
  write_header();
}

/// If desired, this can be called before the first simulation step (but
/// after construction), so that step 0 gets observed.
template <typename ITYP>
void Observable<ITYP>::observe_first()
{
  if (first_observed) return;
  observe();
  first_observed=true;
}

/******************************************************************************
 *
 * SBObservable
 *
 */

/** \class SBObservable
    \brief Step-based observable
    \ingroup Observable
*/
class SBObservable : public Observable<int> {
public:
  SBObservable(SimEnvironment& senv_) : Observable(senv_),
					senv(senv_) {}
protected:
  void step_local();

private:
  SimEnvironment &senv;
} ;

/******************************************************************************
 *
 * SBObservable
 *
 */

/** \class KMCObservable
    \brief Real-time-based observable (as for Kinetic MC)
    \ingroup Observable
*/
class KMCObservable : public Observable<double> {
public:
  KMCObservable(SimEnvironment& env_);

protected:
  void   step_local();
  void   init_local();  
  void   warm_init_local() {Observable<double>::warm_init_local();}

  double obs_time;

private:
  SimEnvironment& env;

  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  virtual void vserial(iarchive_t &ar) {ar >> *this;}
  virtual void vserial(oarchive_t &ar) {ar << *this;}
} ;

inline KMCObservable::KMCObservable(SimEnvironment& env_) :
  Observable(env_),
  obs_time(0.),
  env(env_)
{}

template <typename Archive>
void KMCObservable::serialize(Archive &ar,const unsigned int version)
{
  ar & obs_time;
}


} /* namespace */


BOOST_CLASS_VERSION(glsim::SBObservable,0);
BOOST_CLASS_VERSION(glsim::KMCObservable,0);

#endif /* OBSERVABLE_HH */
