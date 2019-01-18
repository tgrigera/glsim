/*
 * gs_gr.cc -- Compute radial distribution function from a given trayectory
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

/** \file gs_gr.cc
    \ingroup Structure

    \brief Compute the radial distribution function

This is a command-line utility to compute the radial distribution
function from a given trayectory.  It uses the glsim::gr class.
*/

#include "parameters.hh"
#include "gr.hh"
#include "olconfiguration.hh"

struct optlst {
public:
  std::vector<std::string> ifiles;
  int    Nbins;
  double dr;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("gs_gr")
{
  hidden_command_line_options().add_options()
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),
     "input files")
    ;
  command_line_options().add_options()
    ("Nbins,N",po::value<int>(&options.Nbins)->default_value(-1),
     "Number of bins to use (excludes -r)")
    ("deltar,r",po::value<double>(&options.dr)->default_value(-1),
     "With of bins to use (excludes -N)")
    ;

  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options]  ifile [ifile ....]\n\n"
    << "Computes the radial distribution function for the given trajectories,\n"
    << "using the specified bin width or total number of bins.  You must give\n"
    << "either -N or -r.\n"
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

  if (options.Nbins>0 && options.dr>0)
    throw glsim::Runtime_error("cannot give both -N and -r");
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  ifs.read();
  glsim::gr *rdf;
  if (options.Nbins>0)
    rdf=new glsim::gr(conf,options.Nbins);
  else if (options.dr>0)
    rdf=new glsim::gr(conf,options.dr);
  else
    throw glsim::Runtime_error("need one of -N or -r");

  do {
    rdf->push(conf);
  } while (ifs.read());
  std::cout << (*rdf);
  delete rdf;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
