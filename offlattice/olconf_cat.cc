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

    This utility copies one or more configuration or trajectory files
    to a new trajectory file, combining the input files in a single
    output.  It can also copy a subset of the input records.
*/

#include "parameters.hh"
#include "olconfiguration.hh"

static struct {
  std::vector<std::string> ifiles;
  std::string              ofile;
  bool                     id_frame,type_frame;
  long                     first,last,deltarec;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("GS_olconf_cat")
{
  hidden_command_line_options().add_options()
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("ofile,o",po::value<std::string>(&options.ofile)->required(),"Output trajectory file")
    ("first,f",po::value<long>(&options.first)->default_value(0),"First record to copy, negative counts backwards from end")
    ("last,l",po::value<long>(&options.last)->default_value(-1),"Last record to copy, negative counts backwards from end")
    ("increment,i",po::value<long>(&options.deltarec)->default_value(1),
     "Record increment (must be positive, 1 copies all records in range")
    ("id-frame",po::bool_switch(&options.id_frame)->default_value(false),"Assume ids change with time and record as frame variable")
    ("type-frame",po::bool_switch(&options.type_frame)->default_value(false),"Assume types change with time and record as frame variable")
    ;
    positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << " [options] ifile [...]\n\n"
    << "     Copy configurations from configuration or trajectory ifiles to a single\n"
    << "trajectory file.  If you give more than one input file, they will be\n"
    << "considered as one file, respecting the order in which they are named.\n"
    << "You can ask for a subrange of records to be copied (see options -f and -l).\n"
    << "The first record is number 0 and numbering continues consecutively\n"
    << "across file borders.  Negative record numbers count from the end of the file,\n"
    << "with -1 meaning the last record.  Coordinates, velocities and\n"
    << "accelerations will be copied as frame variables if present.\n"
    << "Particle id and type will be copied as header info unless requested\n"
    << "otherwise.\n"
    << "\nOptions:\n";
  show_command_line_options(std::cerr);
  std::cerr  << "\n";
  show_parameters(std::cerr);
  std::cerr  << "\n";
}

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  if (options.first<0) options.first+=ifs.size();
  if (options.last<0) options.last+=ifs.size();

  if (options.first<0) options.first=0;
  if (options.first>ifs.size()) options.first=ifs.size()-1;
  if (options.last<0) options.last=0;
  if (options.last>ifs.size()) options.last=ifs.size()-1;
  
  // Figure out what fields to copy and open output file
  ifs.seek(options.first);
  ifs.read();
  glsim::OLconfiguration        oconf;
  glsim::OLconfig_file::options fopt;
  fopt.time_frame().box_frame();
  if (conf.r) fopt.r_frame();
  if (conf.v) fopt.v_frame();
  if (conf.a) fopt.a_frame();
  if (options.id_frame) fopt.id_frame();
  if (options.type_frame) fopt.type_frame();
  glsim::OLconfig_file of(&conf,fopt);
  of.create(options.ofile.c_str());
  of.write_header();
  
  // cat!
  ifs.seek(options.first);
  while (ifs.pos()<=options.last) {
    ifs.read();
    of.append_record();
    ifs.seek(ifs.pos()+options.deltarec-1);
  } 

}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
