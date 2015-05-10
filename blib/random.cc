/*
 * random.cc -- definitions for random number classes
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

#include <iostream>

#include "random.hh"

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
