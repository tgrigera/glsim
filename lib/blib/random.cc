/*
 * random.cc -- definitions for random number classes
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

#include <iostream>

#include "random.hh"
#include "cerrors.h"

namespace glsim {

// Creation and destruction simply involve locating the corresponding
// scope and eventually calling the corresponding [[gsl_rng]] functions.
// Copying is complicated and not immediatly useful, not supported at
// this time.

Random_number_generator::
Random_number_generator(rng_type type_,const unsigned long seed,
			const char* scope_) :
  scope(scope_),
  type(type_)
{
  if (generator_map.count(scope)>0)
  throw glsim::Logic_error("Only one generator per scope is allowed");

  const gsl_rng_type *gt;
  switch (type) {
  case gsl_rng_mt19937:
    gt=::gsl_rng_mt19937;
    break;
  case gsl_rng_mrg:
    gt=::gsl_rng_mrg;
    break;
  }

  generator_map[scope]=gsl_rng_alloc(gt);
  generator=generator_map[scope];
  glsim_generator_map[scope]=this;
  set_seed(seed);
}

Random_number_generator::Random_number_generator(Random_number_generator& r)
{
  throw glsim::Unimplemented("copying of random number generator not allowed");
}

Random_number_generator& 
Random_number_generator::operator=(Random_number_generator&)
{
  throw glsim::Unimplemented("copying of random number generator not allowed");
}

Random_number_generator::~Random_number_generator()
{
  gsl_rng_free(generator);
  generator_map.erase(scope);
  glsim_generator_map.erase(scope);
}

// Creation and destruction
  
void Random_number_generator::save(std::ostream& os)
{
  void *state=gsl_rng_state(generator);
  os.write((char*) state,gsl_rng_size(generator));
}

void Random_number_generator::load(std::istream& is)
{
  void *state=gsl_rng_state(generator);
  is.read((char*) state,gsl_rng_size(generator));
}

/*
 * Static data
 */

const char* Random_number_generator::default_scope="[default]";
Random_number_generator::genmap_t Random_number_generator::generator_map;
Random_number_generator::glsim_genmap_t Random_number_generator::glsim_generator_map;

}
