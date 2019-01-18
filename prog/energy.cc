/*
 * energy.cc -- generic tool to compute energy of given configurations
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
#include "log.hh"
#include "parameters.hh"
#include "random.hh"
#include "cerrors.h"
#include "olconfiguration.hh"
#include "interactions.hh"

using namespace glsim;

/** \file energy.cc

Generic command-line tool to compute the energy of a set of
configurations.  Needs to be compiled together with another module
that creates the appripriate interaction object.  See ljenergy.cc for
an example.

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
  void show_usage() const;
} ;

CLoptions::CLoptions() : UtilityCL("gs_ljenergy")
{
  hidden_command_line_options().add_options()
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;

  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "\nusage: " << progname << " [options] ifile [..]\n\n"
    << "Computes the energy of the given configurations.\n"
    << "Input files can be trajectory or configuration files.\n"
    << "\nOptions:\n";
  show_command_line_options(std::cerr);
  std::cerr  << "\n";
  show_parameters(std::cerr);
  std::cerr  << "\n";
}

/*****************************************************************************
 *
 * main
 *
 */

Interactions* Interactions_object(OLconfiguration&);

void wmain(int argc,char *argv[])
{
  Interactions *inter;
  CLoptions    opt;
  opt.parse_command_line(argc,argv);

  OLconfiguration conf;
  OLconfig_file cfile(&conf);
  H5_multi_file ifs(options.ifiles,cfile);

  glsim::logs.set_stream(std::cout,glsim::error);
  ifs.read();
  inter=Interactions_object(conf);

  printf("#                  |---- E n e r g y   p e r   p a r t i c l e --|\n");
  printf("#  Step       Time |     Potential         Kinetic          Total|\n");
  double ekin,epot;
  do {
    double P[3];
    inter->fold_coordinates(conf);
    epot=inter->potential_energy(conf)/conf.N;
    ekin=conf.v ? inter->kinetic_energy_and_momentum(conf,P)/conf.N : 0;
    printf("%7ld %10g %15.8e %15.8e %15.8e\n",conf.step,conf.time,epot,ekin,ekin+epot);
  } while (ifs.read());

  delete inter;
}
