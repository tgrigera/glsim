/*
 * Fk.cc -- Compute the intermediate scattering function from a given
 *          trayectory
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

/** \file Fk.cc
    \brief Compute the intermediate scattering function
    \ingroup Structure

This is a command-line utility to compute the intermediate scattering
function (full or only self part).  It uses the scattering function
classes Fk and Fsk.
*/

#include "parameters.hh"
#include "iscatt.hh"
#include "olconfiguration.hh"

struct optlst {
public:
  double k;
  std::vector<std::string> ifiles;
  bool   self_part;
  int    nave;
  long   seed;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage();
} ;

CLoptions::CLoptions() : glsim::UtilityCL("GS_olconf_dump")
{
  command_line_options().add_options()
    ("k",po::value<double>(&options.k)->required(),"wavevector")
    ("self,s",po::bool_switch(&options.self_part)->default_value(false),"compute self part")
    ("nave,n",po::value<int>(&options.nave)->default_value(10),"number of averages")
    ("seed,S",po::value<long>(&options.seed)->default_value(1L),"seed")
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  positional_options().add("k",1).add("ifiles",-1);
}

void CLoptions::show_usage()
{
  std::cerr
    << "usage: " << progname << "[options] k ifile [ifile ....]\n\n"
    << "Computes the (self part of) the intermediate scattering function,"
    << "at the given wavevector and from the given trajectory files.\n"
    << "\n"
    << " Options:\n"
    << "   --self,-s    compute only the self part F_s(k) (default is the full F(k,t)\n"
    << "   --nave,-n nn do nn averages over random directions of the wavevector\n"
    << "                (does not apply for -s, where average is over particles),\n"
    << "                default 10\n"
    << "   --seed,S n   random number seed (with -s used only to generate one random directio)\n"
    << "   --help,-h    show this help\n"
    << "\n";
}


/*****************************************************************************
 *
 * main and deltat
 *
 */

double get_deltat(glsim::H5_multi_file &ifs,glsim::OLconfiguration &conf)
{
  ifs.read();
  double t0=conf.time;
  ifs.read();
  ifs.rewind();
  return conf.time-t0;
}

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);
  glsim::Random_number_generator RNG(glsim::gsl_rng_mt19937,options.seed);

  double deltat=get_deltat(ifs,conf);
  if (options.self_part) {
    glsim::Fsk F(options.k,deltat);
    while (ifs.read()) {
      conf.unfold_coordinates(),
      F.push_config(conf.r,conf.N);
    }
    F.compute_Fsk();
    std::cout << F;
  } else {
    glsim::Fk F(options.k,deltat,options.nave);
    while (ifs.read()) {
      conf.unfold_coordinates(),
      F.push_config(conf.r,conf.N);
    }
    F.compute_Fk();
    std::cout << F;
  }
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
