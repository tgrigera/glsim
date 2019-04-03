/*
 * avesd.cc -- Compute average and standard deviation
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

/** \file avesd.cc
    \ingroup Analysis
    \brief Average and standard deviation from files or stdin

This is a command-line utility to compute the average, standard
deviation and variance of a set of numbers.  It can read from standard
input or from one or several files, in ascii or binary form.

*/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>

#include "parameters.hh"

/*****************************************************************************
 *
 * Options
 *
 */

struct opt {
  bool binary,single_precision;
  bool header;
  std::vector<std::string> files;

  opt() :
    binary(false),
    single_precision(false),
    header(true)
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
    ("binary,b",po::bool_switch(&options.binary)->default_value(false),
     "binary input")
    ("single,s",po::bool_switch(&options.single_precision),
     "if binary input, read input as single precision")
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
    << "Writes average, variance and standard deviation to stdout, computed\n"
    << "using West's recurrence formula.  Reads from stdin if called without\n"
    << "files.\n"
    << "\nOptions:\n";
    show_command_line_options(std::cerr);
}


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

  void           open_next();
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
    break;
  }
}

void read_files::open_next()
{
  if (fnames.size()==0) {
    in=&std::cin;
    return;
  }
  filen++;
  if (filen==fnames.size()) return;
  fin.close();
  fin.open(fnames[filen]);
  in=&fin;
}

void run_recurrence(read_files& in,double &ave,double &var,int &N)
{
  double x;
  double Q,R;

  ave=var=0;
  N=0;
  while (in.read(x)) {
    N++;
    Q=x-ave;
    R=Q/N;
    ave+=R;
    var+=Q*R*(N-1);
  }
  var/=(N-1);
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

  if (options.header) {
    if (options.files.size()>0) {
      printf("# Average from files:\n");
      for (auto s: options.files)
	std::cout << "# " << s << '\n';
    } else {
      std::cout << "# Average from standard input\n";
    }
    std::cout << "#\n#       Average        Variance    StdDeviation  Nsamples\n";
  }

  read_files in(options.files,options.binary,options.single_precision);
  int        nsamples;
  double     ave,var;  

  run_recurrence(in,ave,var,nsamples);
  std::cout << std::setprecision(10) << std::setw(15) << ave
	    << ' ' << std::setw(15) << var
	    << ' ' << std::setw(15) << sqrt(var)
	    << ' ' << std::setw(9) << nsamples << '\n';
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
