/*
 * stochastic.cc --  base environment for stochastic simulations
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

#include "stochastic.hh"

namespace glsim {

void StochasticEnvironment::init_local()
{
  Environment::init_local();

  seed=par.value("stochastic.seed").as<unsigned long>();
  
  delete RNG;
  std::string rng_name=par.value("stochastic.random_number_generator").as<std::string>();
  if (rng_name=="gsl_rng_mt19937")
    RNG=new Random_number_generator(gsl_rng_mt19937,seed,scope_name.c_str());
  else if (rng_name=="gsl_rng_mrg")
    RNG=new Random_number_generator(gsl_rng_mrg,seed,scope_name.c_str());
  else
    throw glsim::Invalid_value(rng_name,"random_number_generator");
}


} /* namespace */
