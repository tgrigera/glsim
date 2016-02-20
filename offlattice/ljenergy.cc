/*
 * ljenergy.cc -- compute energy of given configurations, with LJ potential
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

#include "../config.h"
#include "parameters.hh"
#include "random.hh"
#include "cerrors.h"
#include "olconfiguration.hh"
#include "interactions.hh"
#include "lj.hh"

using namespace glsim;

/** \file ljenergy.cc

Compute energy of a set of configurations, with the LJ potential.

*/

/*****************************************************************************
 *
 * options
 *
 */

static struct ooptions {
  std::vector<std::string>   ifiles;

} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage();
} ;

CLoptions::CLoptions() : UtilityCL("gs_ljenergy")
{
  command_line_options().add_options()
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;

  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage()
{
  std::cerr
    << "\nusage: " << progname << " [options] ifile [..]\n\n"
    << "Computes the energy of the given configurations.\n"
    << "Input files can be trajectory or configuration files.\n"
    << "\nOptions:\n";
  show_parameters(std::cerr);
  show_base_utility_parameters(std::cerr);
  std::cerr  << "\n";
}

/*****************************************************************************
 *
 * main
 *
 */

void wmain(int argc,char *argv[])
{
  CLoptions    opt;
  LennardJones LJ;
  opt.parse_command_line(argc,argv);

  OLconfiguration conf;
  OLconfig_file cfile(&conf);
  H5_multi_file ifs(options.ifiles,cfile);

  while (ifs.read()) {
    Interactions_isotropic_pairwise<LennardJones> inter(LJ,conf);
    printf("%7ld %g %g %d\n",conf.step,conf.time,inter.potential_energy(conf), conf.N);
  }
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}

