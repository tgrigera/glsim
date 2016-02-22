/*
 * gav.cc -- Average in geometrically growing windows
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

#include "parameters.hh"
#include "mfile.hh"
#include "geoave.hh"

/*****************************************************************************
 *
 * Options
 *
 */

struct {
  std::vector<std::string> files;
  bool        header;
  double      base;
  double      t0;
  double      wfactor;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : UtilityCL("gs_gav")
{
  command_line_options().add_options()
    ("base,b",po::value<double>(&options.base)->default_value(1.),
     "set base window length")
    ("t0,t",po::value<double>(&options.t0)->default_value(0.),
     "set time origin")
    ("wfactor,w",po::value<double>(&options.wfactor)->default_value(1.5),
     "set wfactor")
    ;
  hidden_command_line_options().add_options()
    ("input-file",po::value<std::vector<std::string> >(&options.files),
     "input file")
    ;
  positional_options().add("input-file",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "\nusage: " << progname << " [options] [file ... ]\n\n"
    << "     This program takes two-column ascii files (typically time series\n"
    << "files) and averages together all points whose x-coordinate\n"
    << "(we assume it is a time) falls within a window.  The first window starts\n"
    << "at t_0 and is of length base, successive windows grow geometrically\n"
    << "by factor wfactor.  In other words, window boundaries are located at\n\n"
    << "     t_n = t_0 + base * (wfactor^n - 1)/(wfactor-1)\n\n"
    << "     In the output, the time assigned to the average value is the center\n"
    << "of the window, i.e. the abscissas s_n are\n\n"
    << "     s_n = t_n + 0.5 * base * wfactor^n\n\n"
    << "     Input data points may be in any order. The case wfactor=1 is\n"
    << "supported and handled as a special case, yielding windows of fixed\n"
    << "width equal to base. wfactor<1 is not recommended.\n\n"
    << "If no files are given, stdin is used.\n"
    << "\nOptions:\n";
  show_command_line_options(std::cerr);
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
  options.header=!opt.value("terse").as<bool>();
  glsim::Geoave gav(options.t0,options.wfactor,options.base);

  glsim::MFILE mf(options.files);
  char buf[201];
  double time,e;
  while ( !mf.eof() ) {
    fgets(buf,200,mf);
    if (*buf=='#') continue;
    if (sscanf(buf,"%lg %lg",&time,&e)!=2)
      throw glsim::Clib_error(HERE);
    gav.push(time,e);
  }

  if (options.header) {
    if (options.files.empty())
      std::cout << "# Average from standard input\n";
    else {
      std::cout << "# Average from files:\n";
      for (auto f : options.files)
	std::cout << "# " << f << '\n';
    }
    std::cout << "#\n";
  }
  std::cout << gav;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
