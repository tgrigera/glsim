/*
 * random.hh -- provides random number classes
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

#ifndef RANDOM_HH
#define RANDOM_HH


#include <fstream>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/serialization/version.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/base_object.hpp>

#include "exception.hh"

namespace glsim {
  class Random_number_generator;
}

namespace boost { namespace serialization {
template<class Archive>
inline void save_construct_data(Archive & ar,
				const glsim::Random_number_generator *rng,
				const unsigned int file_version);
}}

namespace glsim {

/**  \class Random_number_generator
     \ingroup Random
     \brief The basic random number object

A GSL random number generator (RNG) is created by creating an instance
of this class.  Normally only one `Random_number_generator` object
will exist, which can be shared among many random distributions.  The
internal state of a generator can be saved to a file and later
restored, with the simple interface or through `Boost::serialize`.  We
give here methods to access the raw generator output, but normally the
random numbers will be obtained through one of the distributions
defined below.

*/

/** \ingroup Random
 The different RNGs included in the GSL
*/
enum rng_type {gsl_rng_mt19937,gsl_rng_mrg} ; 

class Random_number_generator {
public:
  Random_number_generator(rng_type type,const unsigned long seed=0,
			  const char* scope=default_scope);
  Random_number_generator(Random_number_generator&);
  Random_number_generator& operator=(Random_number_generator&);
  ~Random_number_generator();
  void set_seed(const unsigned long seed);
  int save(FILE *f);
  int load(FILE *f);
  void save(std::ostream&);
  void load(std::istream&);

  unsigned long raw();
  unsigned long min() const;
  unsigned long max() const;
  unsigned long range() const;

  static const char* default_scope;

private:
  typedef std::map<std::string,gsl_rng*> genmap_t;
  typedef std::map<std::string,Random_number_generator*> glsim_genmap_t;
  static genmap_t generator_map;
  static glsim_genmap_t glsim_generator_map;

  std::string scope;
  rng_type    type;
  gsl_rng     *generator;

  template<typename T> friend class Random_distribution_base;

  friend class boost::serialization::access;
  template <typename Ar>
  void serialize(Ar& ar,const unsigned int version) const;

  template <class Archive>
  friend void boost::serialization::save_construct_data(Archive &,const glsim::Random_number_generator*,
				  const unsigned int);
}  ;


/*
 * The [[friend rdbase]] declarations above refer to
 * [[Random_distribution_base]] specializations that are defined below
 * (outside the namespace).  The last friend is needed to give
 * [[Boost::serialization]] acces to the serialization functions.  We
 * finish the declaration by setting the class version number to keep
 * track of changes for serialization.  This declaration needs to be
 * in the global namespace.
 */


/** Setting the seed and saving/loading to file are immediate,
   saving/loading to a stream are implemented by requesting a binary copy
   of the generator state.
*/
inline void Random_number_generator::set_seed(const unsigned long seed)
{
  gsl_rng_set(generator,seed);
}

inline int Random_number_generator::save(FILE *f)
{
  return gsl_rng_fwrite(f,generator);
}

inline int Random_number_generator::load(FILE *f)
{
  return gsl_rng_fread(f,generator);
}


/// returns an integer in the range [min(),max()].
inline unsigned long Random_number_generator::raw()
{
  return gsl_rng_get(generator);
}

inline unsigned long Random_number_generator::min() const
{
  return gsl_rng_min(generator);
}

inline unsigned long Random_number_generator::max() const
{
  return gsl_rng_max(generator);
}

inline unsigned long Random_number_generator::range() const
{
  return gsl_rng_max(generator)-gsl_rng_min(generator);
}

/// Finally, the serialization methods allow alternative interface to
/// state saving, through the Boost::serialization library.
template <typename Ar>
void Random_number_generator::serialize(Ar& ar,const unsigned int ver) const
{
  boost::serialization::binary_object bo(gsl_rng_state(generator),gsl_rng_size(generator));
  ar & bo;
}

/******************************************************************************/

/** \class Random_distribution_base Random distributions
    \ingroup Random   

The following classes provide random numbers distributed according to
different forms, based on the numbers provided by a
[[Random_number_generator]] object.

The following is the base class for all distributions.  The way to
obtain a random number will be through [[operator()()]], which is pure
virtual at this level.  We provide access to the raw generator's
function here for convenience.
*/
template <typename ranT>
class Random_distribution_base {
public:
  Random_distribution_base(const char* scope=Random_number_generator::default_scope);
  virtual ranT operator()()=0;

  unsigned long raw();
  unsigned long min() const;
  unsigned long max() const;
  unsigned long range() const;

protected:
  gsl_rng                 *generator;
  glsim::Random_number_generator *glsim_generator;

private:
  std::string scope;
} ;

typedef Random_distribution_base<double> rdbase_double;
typedef Random_distribution_base<unsigned long> rdbase_ulong;

template <typename ranT> inline Random_distribution_base<ranT>::
Random_distribution_base(const char* scope_) :
  scope(scope_)
{
  if (Random_number_generator::generator_map.count(scope)==0)
    throw glsim::Invalid_operation("Random number distribution called without Random_number_generator defined in scope "+scope);
  generator=Random_number_generator::generator_map[scope];
  glsim_generator=Random_number_generator::glsim_generator_map[scope];
}

template <typename ranT>
unsigned long Random_distribution_base<ranT>::raw()
{
  return glsim_generator->raw();
}

template <typename ranT>
unsigned long Random_distribution_base<ranT>::min() const {
  return glsim_generator->min();
}

template <typename ranT>
unsigned long Random_distribution_base<ranT>::max() const {
  return glsim_generator->max();
}

template <typename ranT>
unsigned long Random_distribution_base<ranT>::range() const {
  return glsim_generator->range();
}


/*****************************************************************************/

/** \class Uniform_integer
    \ingroup Random
    \brief Uniformly distributed integers

This class generates integer random numbers with a uniform
distribution, ([[unsigned long]]) numbers from 0 to a given integer.

The range will be $[0,\ldots,m-1]$, where $m$ is specified as
argument or taken from the default specified on construction.

*/
class Uniform_integer : public rdbase_ulong {
public:
  Uniform_integer(unsigned long default_m=10,const char* scope=
		  Random_number_generator::default_scope);
  unsigned long operator()();
  unsigned long operator()(unsigned long m);

private:
  unsigned long default_m;
} ;

inline
Uniform_integer::Uniform_integer(unsigned long m,const char *scope) :
  Random_distribution_base(scope),
  default_m(m)
{}

inline unsigned long Uniform_integer::operator()()
{
  return gsl_rng_uniform_int(generator,default_m);
}

inline unsigned long Uniform_integer::operator()(unsigned long m)
{
  return gsl_rng_uniform_int(generator,m);
}

/*****************************************************************************/


/** \class Uniform_real
    \ingroup Random
    \brief Uniformly distributed real numbers in arbitrary range

Uniformly-distributed real numbers in the range $[a,b)$ ($b$
excluded).
*/
class Uniform_real : public rdbase_double {
public:
  Uniform_real(double a=0,double b=1,const char *scope=Random_number_generator::default_scope);
  double operator()();

private:
  double a,range;
} ;

inline Uniform_real::Uniform_real(double a,double b,const char* scope) :
  Random_distribution_base(scope),
  a(a)
{
  range=b-a;
}

inline double Uniform_real::operator()()
{
  return a+range*gsl_rng_uniform(generator);
}

/*****************************************************************************/

/** \class Gaussian_distribution
    \ingroup Random

This class uses GSL RNG to return a real ([[double]]) distributed
according to
\f$$
p(x)\, dx = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2) \,dx.
\f$$
*/
class Gaussian_distribution : public rdbase_double {
public:
  Gaussian_distribution(double sigma=1.,double mean=0.,
			const char* scope=Random_number_generator::default_scope);
  double operator()();

private:
  double sigma,mean;
} ;

inline Gaussian_distribution::Gaussian_distribution(double sigma_,double mean_,
						    const char *scope) :
  Random_distribution_base(scope),
  sigma(sigma_),
  mean(mean_)
{}

inline double Gaussian_distribution::operator()()
{
  return mean + gsl_ran_gaussian(generator,sigma);
}

/*****************************************************************************/

/** \class BivariateGaussian_distribution
    \ingroup Random

This class uses GSL RNG to return two real numbers distributed 
according to the bivariate Gassian distribution with zero mean:
\f$$
 p(x,y) dx dy = {1 \over 2 \pi \sigma_x \sigma_y \sqrt{1-\rho^2}}
 \exp (-(x^2/\sigma_x^2 + y^2/\sigma_y^2 - 2 \rho x y/(\sigma_x\sigma_y))/2(1-\rho^2)) dx dy,
\f$$
where the correlathio coefficient \f$\rho\f$ should lien between -1 and 1.
*/
class BivariateGaussian_distribution : public rdbase_double {
public:
  BivariateGaussian_distribution(double sigmax,double sigmay, double rho,
				 const char* scope=Random_number_generator::default_scope);
  void operator()(double &x,double &y);
  double operator()() {throw Invalid_operation("operator() cannot be used on BivariateGaussian");}

private:
  double sigmax,sigmay,rho;
    
} ;

inline BivariateGaussian_distribution::
BivariateGaussian_distribution(double sx,double sy, double r,const char *scope) :
    Random_distribution_base(scope),
    sigmax(sx), sigmay(sy), rho(r)
{}

inline void BivariateGaussian_distribution::operator()(double &x,double &y)
{
  gsl_ran_bivariate_gaussian(generator,sigmax,sigmay,rho,&x,&y);
}


/*****************************************************************************/

/** \class Exponential_distribution
    \ingroup Random

The exponential distribution is defined by
$$
          p(x)\, dx = {1 \over \mu} \exp(-x/\mu)\, dx.
$$
*/
class Exponential_distribution : public Random_distribution_base<double> {
public:
  Exponential_distribution(double mu,
			   const char *scope=Random_number_generator::default_scope);
  double operator()();

private:
  double mu;
} ;

inline
Exponential_distribution::Exponential_distribution(double mu_,const char *scope) :
  Random_distribution_base(scope),
  mu(mu_)
{}

inline double Exponential_distribution::operator()()
{
  return gsl_ran_exponential(generator,mu);
}

/*****************************************************************************/

/** \class Lognormal_distribution
    \ingroup Random

The log-normal distribution is
$$
p(x)\, dx = {1 \over x \sqrt{2 \pi \sigma^2} } 
          \exp\left[-\frac{(\ln x - \zeta)^2}{2 \sigma^2}\right] \, dx.
$$
*/
class Lognormal_distribution : public rdbase_double {
public:
  Lognormal_distribution(double zeta,double sigma,
			 const char *scope=Random_number_generator::default_scope);
  double operator()();

private:
  double zeta,sigma;
} ;

inline Lognormal_distribution::Lognormal_distribution(double zeta_,
						      double sigma_,
						      const char *scope) :
  Random_distribution_base(scope),
  zeta(zeta_),
  sigma(sigma_)
{}

inline double Lognormal_distribution::operator()()
{
  return gsl_ran_lognormal(generator,zeta,sigma);
}


/*****************************************************************************/
/** \class Levy_distribution
    \ingroup Random

This distribution is defined by
$$
p(x) = {1 \over 2 \pi} \int_{-\infty}^{+\infty} dt \exp(-it x - |c t|^\alpha).
$$
*/
class Levy_distribution : public rdbase_double {
public:
  Levy_distribution(double C,double alpha,
		    const char *scope=Random_number_generator::default_scope);
  double operator()();

private:
  double C,alpha;
} ;

inline Levy_distribution::Levy_distribution(double C_,
					    double alpha_,
					    const char *scope) :
  Random_distribution_base(scope),
  C(C_), alpha(alpha_)
{}

inline double Levy_distribution::operator()()
{
  return gsl_ran_levy(generator,C,alpha);
}

/*****************************************************************************/
/** \class Spherical3d_distribution
    \ingroup Random
    \brief Random direction in three dimensions

This returns a random vector with tip on the unit sphere (random
direction in three dimensions).
*/
class Spherical3d_distribution : public rdbase_double {
public:
  Spherical3d_distribution(const char *scope=Random_number_generator::default_scope) :
  Random_distribution_base(scope) {}

  double operator()();
  void operator()(double *r);
} ;


inline double Spherical3d_distribution::operator()()
{
  throw glsim::Invalid_operation("operator() requires 3 arguments in Spherical3d_distribution");
}

inline void Spherical3d_distribution::operator()(double *r)
{
  gsl_ran_dir_3d(generator,r,r+1,r+2);
}

} /* namespace */



BOOST_CLASS_VERSION(glsim::Random_number_generator,0)


namespace boost { namespace serialization {

///The following is to allow serialization through pointers.
template<class Archive>
inline void save_construct_data(Archive & ar,
				const glsim::Random_number_generator *rng,
				const unsigned int file_version)
{
  ar << rng->scope;
  ar << rng->type;
}

template<class Archive>
inline void load_construct_data(Archive & ar,
				glsim::Random_number_generator *rng,
				const unsigned int file_version)
{
  std::string scope;
  glsim::rng_type type;
  ar >> scope;
  ar >> type;
  ::new(rng) glsim::Random_number_generator(type,0,scope.c_str());
}

}} 


#endif /* RANDOM_HH */
