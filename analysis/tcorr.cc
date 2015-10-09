/*
 * tcorr.cc -- Read time series and compute time correlations
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

/** \file  tcorr.cc
    \ingroup TCorr
    \brief Read time series and compute time correlations

This programs use the functions defined in the Analysis module of the
library to compute the time autocorrelation of real or complex
quantities, reading two- or three-column input (\f$t\f$ and \f$A\f$),
from a named file or standard input.  Input/output makes up for most
of the code.

*/

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <errno.h>

#include <iostream>
#include <vector>
#include <complex>

#include "config.h"
#include "parameters.hh"
#include "timecorr.hh"
#include "exception.hh"
#include "mfile.hh"

#ifdef HAVE_LIBFFTW3
typedef glsim::RealFFTW     rFFT;
typedef glsim::ComplexFFTW  cFFT;
#else
typedef glsim::RealFFT_gsl_2n     rFFT;
typedef glsim::ComplexFFT_gsl_2n  cFFT;
#endif /* HAVE_LIBFFTW3 */

/*****************************************************************************
 *
 * Options
 *
 */

struct o {
  std::vector<std::string> ifiles;
  bool   read_times;
  bool   complex_data;
  bool   connected;
  bool   normalize;
  bool   t0set;
  double t0;
  int    fractdiff;
  double deltat;
  bool   use_fft;
  bool   fixed_tw;

  o() :
    read_times(true),
    deltat(0),
    t0(0),
    t0set(false)
  {}

} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage();

private:
  static void handle_deltat(double),handle_t0(double),handle_w(bool);
} ;

CLoptions::CLoptions() : UtilityCL("gs_tcorr")
{
  command_line_options().add_options()
    ("input-file",po::value<std::vector<std::string> >(&options.ifiles),
     "input file")
    ("fft,f",po::bool_switch(&options.use_fft)->default_value(false),
     "use FFT method")
    ("connected,C",po::bool_switch(&options.connected)->default_value(false),
     "compute connected correlation")
    ("normalize,N",po::bool_switch(&options.normalize)->default_value(false),
     "normalize connected correlation")
    ("frac,n",po::value<int>(&options.fractdiff)->default_value(2),
     "report up to npoints/n correlation points")
    ("delta-t,D",po::value<double>(&options.deltat)->notifier(handle_deltat),
     "set time increment")
    ("time-origin,t",po::value<double>(&options.t0)->notifier(handle_t0),
     "set time origin")
    ("fixed-tw,w",po::bool_switch(&options.fixed_tw)->default_value(false)->notifier(handle_w),
     "fixed tw")
    ("complex,x",po::bool_switch(&options.complex_data)->default_value(false),
     "read complex data")
    ;
  positional_options().add("input-file",-1);
}

void CLoptions::handle_deltat(double dum)
{
  if (options.t0set)
    throw glsim::Runtime_error("cannot give both -D and -t");
  options.read_times=false;
}

void CLoptions::handle_t0(double dum)
{
  if (!options.read_times)
    throw glsim::Runtime_error("cannot give both -D and -t");
  options.t0set=true;
}

void CLoptions::handle_w(bool fixed)
{
  options.fixed_tw=fixed;
  if (fixed && options.use_fft)
    throw glsim::Runtime_error("cannot give -w with FFT method (-f)");
}

void CLoptions::show_usage()
{
  std::cerr 
    << "\nusage: " << progname << "[options] [file ...]\n\n"
    << "Computes the time-autocorrelation from the two- or three-column files or\n"
    << "from stdin and writes to stdout.  Unless -w is given, TTI is assumed\n"
    << "(i.e. will average over all time origins.\n\n"
    << " Options:\n"
    << "   -h, --help            Show this help\n"
    << "   -f, --fft             Use the FFT method (default is the direct O(N^2) method)\n"
    << "   -C, --connected       Compute connected correlations\n"
    << "   -N, --normalize       Normalize to C(0)=1 (only good with -C)\n"
    << "   -n arg, --frac        Compute only (number of data)/arg time differences (must be\n"
    << "                         an integer, default arg=2)\n"
    << "   -D arg, --delta-t     Set time step (Delta t) to arg and read only ONE column (two with -x)\n"
    << "   -t arg, --time-origin Set time origin (discards all data with time<arg)\n"
    << "   -w, --fixed-tw        Compute corr at FIXED tw (equal to time given with -t)\n"
    << "   -x, --complex         Observable is COMPLEX; expects three columns (two with -D)\n"
    << "\n";

  // std::cerr << command_line_options();

}

/*****************************************************************************
 *
 * Read data
 *
 */

void read_data(glsim::MFILE &f,std::vector<double>& a)
{
  char   buf[200];
  double t,at;

  if (options.read_times) {
    int n=0;
    while ( !f.eof() ) {
      fgets(buf,200,f);
      if (*buf=='#') continue;
      if (sscanf(buf,"%lg %lg",&t,&at)!=2)
          throw glsim::Clib_error(HERE);
      if (options.t0set && t<options.t0) continue;
      a.push_back(at);
      if (n==0) options.deltat=t;
      if (n==1) options.deltat=t-options.deltat;
      n++;
    }
  } else {
    while ( !f.eof() ) {
      fgets(buf,200,f);
      if (*buf=='#') continue;
      sscanf(buf,"%lg",&at);
      if (errno) throw glsim::Clib_error(HERE);
      a.push_back(at);
    }
  }
}

void read_data(glsim::MFILE &f,glsim::vcomplex& a)
{
  char   buf[200];
  double t,atr,ati;

  if (options.read_times) {
    int n=0;
    while ( !f.eof() ) {
      fgets(buf,200,f);
      if (*buf=='#') continue;
      sscanf(buf,"%lg %lg %lg",&t,&atr,&ati);
      if (errno) throw glsim::Clib_error(HERE);
      if (options.t0set && t<options.t0) continue;
      a.push_back(glsim::dcomplex(atr,ati));
      if (n==0) options.deltat=t;
      if (n==1) options.deltat=t-options.deltat;
      n++;
    }
  } else {
    while ( !f.eof() ) {
      fgets(buf,200,f);
      if (*buf=='#') continue;
      sscanf(buf,"%lg %lg",&atr,&ati);
      if (errno) throw glsim::Clib_error(HERE);
      a.push_back(glsim::dcomplex(atr,ati));
    }
  }
}

/*****************************************************************************
 *
 * main
 *
 */

void wmain(int argc,char *argv[])
{
  rFFT                   *realft=0;
  cFFT                   *complexft=0;
  std::vector<double>    a,corr;
  const std::vector<double> *acorr=0;
  glsim::vcomplex        ac,corrc;
  const glsim::vcomplex *acorrc=0;
  double                deltat;

  CLoptions opt;
  opt.parse_command_line(argc,argv);

  glsim::MFILE fin(options.ifiles);
  if (options.complex_data)
    read_data(fin,ac);
  else
    read_data(fin,a);
  int nt=a.size()/options.fractdiff;
  if (nt<2) throw glsim::Runtime_error("Too few points!");

  /* Call function to compute requested correlation */
  if (options.use_fft) {
    if (options.complex_data) {
      complexft=new cFFT(glsim::FFT::in_place);
      if (options.connected)
	glsim::correlation_connected_1d_tti_fft(ac,*complexft,nt,options.normalize);
      else
	glsim::correlation_1d_tti_fft(ac,*complexft,nt);
      acorrc=&complexft->tdata();
    } else {
      realft=new rFFT(glsim::FFT::in_place);
      if (options.connected)
  	correlation_connected_1d_tti_fft(a,*realft,nt,options.normalize);
      else
  	correlation_1d_tti_fft(a,*realft,nt);
      acorr=&realft->tdata();
    }
  } else {
    if (options.complex_data) {
      if (options.connected)
	glsim::correlation_connected_1d_tti_direct(ac,corrc,nt,options.normalize);
      else
	glsim::correlation_1d_tti_direct(ac,corrc,nt);
      acorrc=&corrc;
    } else {
      if (options.fixed_tw) 
  	if (options.connected)
	  glsim::correlation_connected_1d_tw_direct(a,corr,0,options.normalize);
  	else
	  glsim::correlation_1d_tw_direct(a,corr,0);
      else
  	if (options.connected)
	  glsim::correlation_connected_1d_tti_direct(a,corr,nt,options.normalize);
  	else
	  glsim::correlation_1d_tti_direct(a,corr,nt);
      acorr=&corr;
    }
  }

  /* Integral */
  double sum=0;
  for (int i=0; i<acorr->size(); i++)
    sum+=(*acorr)[i];
    sum*=options.deltat;

  /* Output */
  printf("# Time correlation\n#\n");
  printf("#\n# Integral = %g\n#\n",sum);
  if (options.connected)
    printf("# time    <[a(0)-<a>][a(t)-<a>]>\n");
  else
    printf("# time    <a(0)a(t)>\n");
  if (options.complex_data)
    for (int i=0; i<nt; i++)
      printf("%g %g %g\n",i*options.deltat,(*acorrc)[i].real(),(*acorrc)[i].imag());
  else
    for (int i=0; i<nt; i++)
      printf("%g %g\n",i*options.deltat,(*acorr)[i]);

  delete realft,complexft;
}


int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
