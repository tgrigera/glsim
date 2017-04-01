/*
 * olconf_dump.cc -- dump configurations summary data to stdout (ascii)
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

/** \file olconf_dump.cc
    \ingroup Offlattice
    \brief Dump summary data on given configurations to stdout

    This is so far somewhat primitive.  However, for full dump you can
    use H5_dump.
*/

#include "parameters.hh"
#include "olconfiguration.hh"

static struct {
  bool   dump_positions;
  std::vector<std::string> ifiles;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("gs_olconf_dump")
{
  hidden_command_line_options().add_options()
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("positions,p",po::bool_switch(&options.dump_positions),
     "Dump positions in ASCII")
    ;
  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "ifile [ifile ....]\n\n"
    << "Dump information about configurations\n"
    << "\n";
  show_command_line_options(std::cerr);
}

void dump_positions(glsim::H5_multi_file &ifs,glsim::OLconfiguration &conf)
{
  printf("#  Step   Id            x            y            z\n");
  while (ifs.read()) {
    for (int i=0; i<conf.N; ++i)
      printf("%7ld %4d %12.6e %12.6e %12.6e\n",conf.step,conf.id[i],
	     conf.r[i][0],conf.r[i][1],conf.r[i][2]);
  }
}

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file ifs(options.ifiles,cfile);

  if (options.dump_positions) dump_positions(ifs,conf);
  else {
    printf("# Step  Time  N  density\n");
    while (ifs.read()) {
      double vol=conf.box_length[0]*conf.box_length[1]*conf.box_length[2];
      printf("%7ld %g  %d  %g\n",conf.step,conf.time,conf.N,(double) conf.N/vol);
    }
  }
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
