/*
 * Fk.cc -- Compute the intermediate scattering function from a given
 *          trayectory
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
 * help convince funding agencies that it is worth funding.
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
  bool   substract_cm;
  int    nave;
  long   seed;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("GS_olconf_dump")
{
  hidden_command_line_options().add_options()
    ("k",po::value<double>(&options.k)->required(),"wavevector")
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("substract-cm,c",po::bool_switch(&options.substract_cm)->default_value(false),
     "substract center of mass at each time (assumes all particles have the same mass)")
    ("self,s",po::bool_switch(&options.self_part)->default_value(false),
     "compute only the self part F_s(k) [default is the full F(k,t)]")
    ("nave,n",po::value<int>(&options.nave)->default_value(10),
     "do arg averages over random directions of the wavevector")
    ("seed,S",po::value<long>(&options.seed)->default_value(1L),
     "random number seed (with -s used to generate only one random direction)")
    ;

  positional_options().add("k",1).add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] k ifile [ifile ....]\n\n"
    << "Computes the (self part of) the intermediate scattering function,"
    << "at the given wavevector and from the given trajectory files.\n"
    << "\n"
    << " Options:\n";
  show_command_line_options(std::cerr);
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

void substract_cm(glsim::OLconfiguration &conf)
{
  std::vector<double> CM=conf.center_of_mass();
  for (int i=0; i<conf.N; ++i) {
    conf.r[i][0]-=CM[0];
    conf.r[i][1]-=CM[1];
    conf.r[i][2]-=CM[2];
  }
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
      conf.unfold_coordinates();
      if (options.substract_cm) substract_cm(conf);
      F.push_config(conf.r,conf.N);
    }
    F.compute_Fsk();
    std::cout << F;
  } else {
    glsim::Fk F(options.k,deltat,options.nave);
    while (ifs.read()) {
      conf.unfold_coordinates();
      if (options.substract_cm) substract_cm(conf);
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
