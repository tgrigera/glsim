/*
 * histocl.cc -- Command-line utility to create histograms
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
 * help convince funding agencies that it is worth funding.  * glsim
 * distribution.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */


/** \file histocl.cc
    \ingroup Analysis
    \brief Histogram from files or stdin

This is a command-line utility to create histograms.  It can read from
standard input or from one or several files, in ascii or binary form.

*/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>

#include "parameters.hh"
#include "histogram.hh"

/*****************************************************************************
 *
 * Options
 *
 */

struct opt {
  int         nbins;
  std::string smin,smax;
  double      min,max;
  bool        probability;
  bool        binary,single_precision;
  std::vector<std::string> files;

  opt() :
    binary(false),
    single_precision(false),
    probability(false)
  {}
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : UtilityCL("gs_avesd")
{
  command_line_options().add_options()
    ("nbins",po::value<int>(&options.nbins)->required())
    ("min",po::value<std::string>(&options.smin)->required())
    ("max",po::value<std::string>(&options.smax)->required())
    ("input-file",po::value<std::vector<std::string> >(&options.files),
     "input file")
    ("binary,b",po::bool_switch(&options.binary)->default_value(false),
     "binary input")
    ("probability,p",po::bool_switch(&options.probability)->default_value(false),
     "probability output")
    ("single,s",po::bool_switch(&options.single_precision),
     "read binary input as single precision")
    ;
  positional_options().add("nbins",1).add("min",1).add("max",1).add("input-file",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "\nusage: " << progname << " [options] nbins min max [file ... ]\n\n"
    << "Computes a histogram from data given in files (or stdin if called without\n"
    << "files.  Writes to stdout.  Data outside the range [min,max] will be counted\n"
    << "separately as outliers.\n"
    << "If you request probability output (-p), be aware that outliers are taken into\n"
    << "account when computing total area, i.e. total area is 1 if there are no outliers,\n"
    << "otherwise it is equal to the fraction of non-outliers.\n"
    << "\nOptions:\n"
    << "   -b, --binary Binary input\n"
    << "   -s, --single If binary input, read as single precision (float)\n"
    << "                otherwise assume double precision\n"
    << "   -h, --help   Show this help\n"
    << "   -p, --probabily Probabilty output (i.e. normalize so that total area is 1).\n"
    << "\nNegative numbers need to be quoted (as \"-10\").  Depending on your shell, you\n"
    << "will need to escape the quotes so that they are passed to the program, in bash\n"
    << "you must write \\\"-10\\\".\n"
    << '\n';
}

/*******************************************************************************
 *
 * read_files
 *
 */

class read_files {
public:
  read_files(std::vector<std::string> &files,bool binary_,bool single_precision_);
  bool read(double&);

private:
  bool                     binary,single_precision;
  std::vector<std::string> &fnames;
  int                      filen;
  std::istream             *in;
  std::ifstream            fin;

  bool          open_next();
} ;

read_files::read_files(std::vector<std::string> &files,bool binary_,
		       bool single_precision_) :
  fnames(files),
  binary(binary_),
  single_precision(single_precision_),
  filen(-1),
  in(0)
{
  open_next();
}

bool read_files::read(double &x)
{
  switch (binary) {
  case false:
    *in >> x;
    if (!in->eof()) return true;
    open_next();
    *in >> x;
    return !in->eof();
    break;
  case true:
    in->get();
    if (in->eof()) open_next();
    else in->unget();

    if (single_precision) {
      float d;
      in->read((char *) &d,sizeof(d));
      x=d;
    } else
      in->read((char *) &x,sizeof(x));
    return !in->eof();
  }
}

bool read_files::open_next()
{
  if (fnames.size()==0) {
    in=&std::cin;
    return true;
  }
  filen++;
  if (filen==fnames.size()) return false;
  fin.close();
  fin.open(fnames[filen]);
  in=&fin;
  return true;
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

  if (options.files.size()>0) {
    std::cout << "# Histogram from files:\n";
    for (auto s: options.files)
      std::cout << "# " << s << '\n';
  } else {
    std::cout << "# Histogram from standard input\n";
  }
  if (options.smin[0]=='\"') {
    options.smin.erase(0,1);
    options.smin.erase(options.smin.size()-1,1);
  }
  options.min=std::stod(options.smin);
  if (options.smax[0]=='\"') {
    options.smax.erase(0,1);
    options.smax.erase(options.smin.size()-1,1);
  }
  options.max=std::stod(options.smax);

  read_files in(options.files,options.binary,options.single_precision);
  glsim::Histogram histog(options.nbins,options.min,options.max);

  double x;
  while (in.read(x)) histog.push(x);
  
  std::cout << "#\n# Total points processed = " << histog.npoints() << '\n'
	    << "# Outliers total         = " << histog.outliers() << '\n'
	    << "# Outliers below         = " << histog.outliers_low() << '\n'
	    << "# Outliers above         = " << histog.outliers_high() << '\n'
	    << "# Median (coarse-grained)= " << histog.median() << '\n';
  std::cout << "# Bin  Frequency\n";
  histog.probability_output(options.probability);
  std::cout << histog;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
