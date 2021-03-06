/*
 * tau.cc -- Estimate relaxation time from a connected correlation
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

/** \file  tau.cc
    \ingroup TCorr
    \brief Estimate correlation time from a connected correlation function

Compute autocorrelation times given normalized connected correlations,
using Sokal method.
*/

#include <iostream>
#include <cstdio>
#include <vector>
#include <list>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "exception.hh"
#include "parameters.hh"
#include "avevar.hh"

/*****************************************************************************
 *
 * Options
 *
 */

struct o {
public:
  std::vector<std::string> files;
  double alpha;
  enum  {none, tau, alphatest, integraltest} job;
  bool   header;

  o() :
    alpha(-1),
    job(none),
    header(true)
  {}
} options;


class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;

private:
  static void handle_alpha(double),handle_alphatest(double);
  static void handle_integraltest(bool)  ;
} ;

CLoptions::CLoptions() : UtilityCL("gs_tcorr")
{
  command_line_options().add_options()
    ("alpha,a",po::value<double>(&options.alpha)->notifier(handle_alpha),
     "Set aplha.  Output will be tau with jacknife error")
    ("alpha-test,t",po::value<double>(&options.alpha)->notifier(handle_alphatest),
     "Do alpha test with aplpha=1,2,...,arg")
    ("integral-test,i",po::bool_switch()->notifier(handle_integraltest),
     "Do integral test")
    ;
  hidden_command_line_options().add_options()
    ("input-file",po::value<std::vector<std::string> >(&options.files),
     "input files");
  positional_options().add("input-file",-1);
}

void CLoptions::handle_alpha(double alpha)
{
  options.job=o::tau;
}

void CLoptions::handle_alphatest(double alpha)
{
  options.job=o::alphatest;
}

void CLoptions::handle_integraltest(bool i)
{
  if (i) options.job=o::integraltest;
  if (options.job==o::none)
    throw glsim::Runtime_error("Must give -i or -a or -t with positive argument\n");
}

void CLoptions::show_usage() const
{
  std::cerr 
    << "\nusage: " << progname << "[options] file [...]\n\n"
    << "Computes the relaxation time with Sokal's method. Input files must\n"
    << "be different samples of a connected, normalized correlation function\n"
    << "(two-column data as written by tcorr)\n\n"
    << "Options:\n";
  show_command_line_options(std::cerr);
}

/*****************************************************************************
 *
 * Compute tau
 *
 */

typedef std::vector<double> corr_t;
typedef std::list<corr_t*> corr_list_t;

double tau_1samp(glsim::AveVar_vector &c,double alpha,double deltat)
{
  int M=0;
  double tau=c.ave(0)/2;
  while (M+.5<alpha*tau && M<c.size()-1) tau+=c.ave(++M);
  if (M>=c.size()-1) std::cerr << "WARNING: M overflow\n";
  return tau*deltat;
}

double tau_nsamp(corr_list_t &c,double alpha,double deltat)
{
  glsim::AveVar_vector cav;
  corr_list_t::iterator i;

  for (i=c.begin(); i!=c.end(); ++i )
    cav.push(**i);

  return tau_1samp(cav,alpha,deltat);
}

void tau_jacknife(corr_list_t &corr,double deltat)
{
  double tauave=tau_nsamp(corr,options.alpha,deltat);
  int n=corr.size();

  if (n==1) { /* can't do jacknife */
    if (options.header) std::cout << "# tau [no error estimate with single sample]\n";
    std::cout << tauave  << '\n';
    return;
  }

  std::vector<double> tau;
  for (int i=0; i<n; i++) {
    corr_t *c=corr.front();
    corr.pop_front();
    tau.push_back(tau_nsamp(corr,options.alpha,deltat));
    corr.push_back(c);
  }

  /* compute jacknife variance estimator */
  double v=0;
  for (auto vi : tau)
    v+=(vi-tauave)*(vi-tauave);
  v*=(tau.size()-1.)/tau.size();

  if (options.header) std::cout << "# tau dtau\n";
  std::cout << tauave << "  " << sqrt(v) << '\n';
}

void tau_alphatest(corr_list_t &corr,double alphamax,double deltat)
{
  std::cout << "# alpha tau\n";
  for (double alpha=1; alpha<=alphamax; alpha+=1.) {
    double tau=tau_nsamp(corr,alpha,deltat);
    std:: cout << alpha << "  " << tau << '\n';
  }
}

void tau_integraltest(corr_list_t &corr,double deltat)
{
  glsim::AveVar_vector cav;
  for (corr_list_t::iterator c=corr.begin(); c!=corr.end(); ++c )
    cav.push(**c);

  std::cout << "# y   <\\int_0^y C(t) dt>  Nsamp\n";
  double in=0.5*deltat*cav.ave(0);
  std::cout << 0.5*deltat << "  " << in << "  " << cav.N(0) << '\n';
  for (int M=1; M<cav.size(); M++) {
    in+=deltat*cav.ave(M);
    std::cout << (M+.5)*deltat << "  " << in << "  " << cav.N(M) << '\n';
  }
}


/*****************************************************************************
 *
 * Read data
 *
 */

void read_data(const char *fname,std::vector<double>& a,double &deltat)
{
  char   buf[200];
  double t,at;

  FILE *f=fopen(fname,"r");
  if (!f) throw glsim::Clib_error(HERE);
  int n=0;
  while ( ungetc(fgetc(f),f)!=EOF ) {
    fgets(buf,200,f);
    if (*buf=='#') continue;
    if (sscanf(buf,"%lg %lg",&t,&at)!=2)
      throw glsim::Clib_error(HERE);
    a.push_back(at);
    if (n==0) deltat=t;
    if (n==1) deltat=t-deltat;
    n++;
  }
  fclose(f);
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

  corr_t      *c;
  corr_list_t corr;
  double      deltat;

  for (auto file : options.files) {
    c=new corr_t;
    read_data(file.c_str(),*c,deltat);
    corr.push_back(c);
  }

  switch (options.job) {
  case o::alphatest:
    tau_alphatest(corr,options.alpha,deltat);
    break;
  case o::integraltest:
    tau_integraltest(corr,deltat);
    break;
  case o::tau:
    tau_jacknife(corr,deltat);
    break;
  }

  while (!corr.empty()) {
    delete corr.front();
    corr.pop_front();
  }
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
