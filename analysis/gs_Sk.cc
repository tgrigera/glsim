/*
 * gs_Sk.cc -- Compute static structure factor from a given trayectory
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

/** \file gs_Sk.cc
    \brief Compute the static structure factor
    \ingroup Structure

This is a command-line utility to compute the static structure factor, either along a given direction or isotropically. It uses the structure factor class Sk.
*/

#include "parameters.hh"
#include "Sk.hh"
#include "olconfiguration.hh"

struct optlst {
public:
  std::vector<std::string> ifiles;
  int  Nk;
  bool isotropic;
  int  kdir;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("gs_Sk")
{
  hidden_command_line_options().add_options()
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("Nk,N",po::value<int>(&options.Nk)->default_value(-1),
     "Number of wavevectors to compute (start is k=0, increment is given by box length).  If negative, try to guess a suitable value from the density.")
    ("isotropic,i",po::bool_switch(&options.isotropic)->default_value(false),
     "compute isotropic S(k) (time is O(N))")
    ("kdirection,k",po::value<int>(&options.kdir)->default_value(0),
     "direction of the k vector (only along the axes, 0=X, 1=Y, 2=Z)")
    ;

  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] ifile [ifile ....]\n\n"
    << "Computes the static structure factor for evenly spaced k vectors.\n"
    << "\n"
    << " Options:\n";
  show_command_line_options(std::cerr);
}


/*****************************************************************************
 *
 * main
 *
 */

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  ifs.read();
  if (options.Nk<0) {
    options.Nk=(int) ceil(conf.box_length[0]*pow(conf.number_density(),1./3.));
    options.Nk*=4;
  }
  glsim::Sk Sk(conf.box_length,options.Nk);
  do {
    if (options.isotropic)
      Sk.push(conf);
    else
      Sk.push(conf,options.kdir);
  } while (ifs.read());
  std::cout << Sk;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
