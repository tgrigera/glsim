/*
 * olconf_cat.cc -- copy configurations
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

/** \file olconf_cat.cc
    \ingroup Offlattice
    \brief A cat program for OLconfiguration files

    This is so far very primitive...
*/

#include "parameters.hh"
#include "olconfiguration.hh"

static struct {
  std::string ifile;
  std::string ofile;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("GS_olconf_cat")
{
  hidden_command_line_options().add_options()
    ("ifile",po::value<std::string>(&options.ifile)->required(),"input file")
    ("ofile",po::value<std::string>(&options.ofile)->required(),"output file")
    ;
  positional_options().add("ifile",1).add("ofile",1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "ifile ofile\n\n"
    << "Copy configurations from ifile to ofile\n"
    << "\n";
}

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration cin;
  cin.load(options.ifile);
  glsim::OLconfiguration cout(cin);
  cout.save(options.ofile);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
