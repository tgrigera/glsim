/*
 * traj2h5md.cc -- convert OLconfiguration trajectories to H5MD format
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
#include "h5md.hh"

using namespace glsim;

/** \file traj2h5md.cc

This provides a utility (GS_traj2h5md) to produce trajectory files in
H5MD format.  So far rather primitve and lacking options.

*/

/*****************************************************************************
 *
 * options
 *
 */

static struct ooptions {
  std::string   traj_file;
  std::string   h5md_file;
  std::string   author;
  std::string   email;

  ooptions() :
    email("") {}
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage();
} ;

CLoptions::CLoptions() : UtilityCL("GS_olconf_create")
{
  command_line_options().add_options()
    ("traj_file",po::value<std::string>(&options.traj_file)->required(),"input trajectory file")
    ("h5md_file,o",po::value<std::string>(&options.h5md_file)->required(),"output file in h5md format")
    ("author",po::value<std::string>(&options.author)->required(),"author name")
    ("email",po::value<std::string>(&options.email),"author email")
    ;
  positional_options().add("traj_file",1);
}

void CLoptions::show_usage()
{
  std::cerr
    << "\nusage: " << progname << " [options] trajfile [..]\n\n"
    << "Reads from trajfiles or configuration files (in glsim OLconfig_file format)\n"
    << "and writes a trajectory file in H5MD format.\n"
    << "So far is very primitive, writing only positions vs time, and does not support\nevolving box size, types etc.\n"
    << "\nOptions:\n"
    << "  -o h5mdfile      Output file [REQUIRED]\n"
    << "  --author string  Author name [REQUIRED]\n"
    << "  --email  string  Author email\n"
    << "\n";
}

/*****************************************************************************
 *
 * main
 *
 */

void wmain(int argc,char *argv[])
{
  CLoptions opt;
  opt.parse_command_line(argc,argv);

  OLconfiguration conf;
  OLconfig_file cfile(options.traj_file.c_str(),&conf);
  cfile.open();
  cfile.read_header();
  cfile.read_record(0);

  H5MD h5md;
  h5md.create(options.h5md_file.c_str(),conf.N,
	      options.author.c_str(),options.email.c_str());
  h5md.set_box(conf.box_length);

  for (int n=0; n<cfile.size(); ++n) {
    cfile.read_record(n);
    h5md.append_positions(conf);
  }
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
