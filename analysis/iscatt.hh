/*
 * iscatt.hh -- Classes to compute the intermediate scattering function
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

#ifndef ISCATT_HH
#define ISCATT_HH


/** \class section{Intermediate scattering function}

Compute the intermediate scattering function $F(k,t)$ and self part of
intermediate scattering function $F_s(k,t)$,
$$ F(k,t) = \frac{1}{N}\langle \rho_k(t)\rho_{-k}(0)\rangle =
\frac{1}{N} \sum_{i,j}\langle e^{i k\cdot[r_i(t)-r_j(0)]},$$
$$ F_s(k,t) = \frac{1}{N} \sum_{i}\langle e^{i k\cdot[r_i(t)-r_i(0)]}
= \langle e^{-i k \cdot r_i(t)} e^{i k\cdot r_i(0)}.$$

\subsection{Through time correlations}

The following functions use the time autocorrelation functions from
the analysis library, and average over $k$ taking a user-given number
of random directions.
*/

#include "config.h"
#include "random.hh"
#include "timecorr.hh"

#ifdef HAVE_LIBFFTW3
typedef glsim::RealFFTW     rFFT;
typedef glsim::ComplexFFTW  cFFT;
#else
typedef glsim::RealFFT_gsl_2n     rFFT;
typedef glsim::ComplexFFT_gsl_2n  cFFT;
#endif /* HAVE_LIBFFTW3 */

namespace glsim {

class Fk {
public:
  Fk(double k,double deltat,int Nav);

  Fk& push_config(double r[][3],int N);
  Fk& compute_Fk();
  const vcomplex& Fk_data() const {return Fk_;} 

private:
  int                               Nav,Npart;
  double                            k,deltat;
  std::vector<std::vector<double> > kr;
  std::vector<vcomplex>             rhok_;
  vcomplex                          Fk_;
  
  dcomplex rho_k(double r[][3],double k[]);

  friend std::ostream& operator<<(std::ostream&,const Fk&); 
} ;


class Fsk {
public:
  Fsk(double k,double deltat);

  Fsk& push_config(double r[][3],int N);
  Fsk& compute_Fsk();
  const vcomplex& Fsk_data() const {return Fsk_;}

private:
  int                   Npart;
  double                k,deltat;
  std::vector<double>   kr;
  std::vector<vcomplex> expkr;
  vcomplex              Fsk_;

  friend std::ostream& operator<<(std::ostream&,const Fsk&); 
} ;


} /* namespace */

#endif /* ISCATT_HH */














#ifdef KLNLKNLKN



// For the self part ($F_s(k)$) we take only one random direction, as
// we average over the particles.

void Fk::Fsk(vcomplex& Fsk)
{
  spherical_random_3d sr;
  double kr[3];
  sr(kr);
  for (int i=0; i<3; i++) kr[i]*=k;

  std::vector<vcomplex> expkr(conf.N);

  while (tf.read(conf))
    for (int i=0; i<conf.N; i++)
      expkr[i].push_back( exp( tcomplex(0,-1)*
			       (kr[0]*conf.r[i][0] + kr[1]*conf.r[i][1] + 
				kr[2]*conf.r[i][2]) ) );

  // compute Fk for each direction and average
  tcomplex fac=tcomplex(1./conf.N,0);
  int Fklen=expkr[0].size()/2;
  Fsk.clear();
  Fsk.resize(Fklen);
  cFFT FF(FFT::in_place);
  for (int i=0; i<conf.N; i++) {
    correlation_1d_tti_fft(expkr[i],FF,Fklen);
    for (int j=0; j<FF.tdata_rw().size(); j++) Fsk[j]+=fac*FF.tdatum(j);
  }
}

@ A full program to compute $F(k,t)$ using the routines above

<<Fk.cc>>=
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <assert.h>

#include "glsim/errdeb.h"
#include "offlattice/olconfig.hh"
#include "offlattice/iscatt.hh"
<<Fk.cc options>>
<<Fk.cc write>>
<<Fk.cc main>>

<<Fk.cc main>>=
int main(int argc,char *argv[])
{
  try {

    options opt(argc,argv);
    random_number_generator RNG(gsl_rng_mt19937,opt.seed);
    traj_istream *tf;
    if (opt.gromacs_trajectory)
      tf=new traj_istream_gromacs(opt.ifiles);
    else
      tf=new traj_istream_netcdf(opt.ifiles);
  
    // determine deltat
    double deltat;
    olconfig *oc;
    oc=new olconfig("");
    assert(tf->read(*oc));  // sometimes the first record has strange time
    assert(tf->read(*oc));
    deltat=oc->time;
    assert(tf->read(*oc));
    deltat=oc->time-deltat;
    delete oc;
    tf->rewind();

    vcomplex Fk;
    switch(opt.which_Fk) {
    case options::full:
      Fk_num_ave(*tf,opt.k,Fk,opt.nave);
      break;
    case options::self:
      Fsk_num_ave(*tf,opt.k,Fk);
      break;
    }
    Fk_write(opt,deltat,Fk);

    delete tf;
  } catch (std::exception &e) {
    std::cerr << "Exception: " << e.what() << '\n';
  } catch (ex::exception &e) {
    std::cerr << e.what() << " (at: " << e.where() << ")\n";
  }

  return 0;
}

<<Fk.cc write>>=
void Fk_write(options &opt,double deltat,vcomplex &Fk)
{
  switch(opt.ave_type) {
  case options::numerical:
    switch(opt.which_Fk) {
    case options::full:
      printf("# time   Fk'   Fk''\n");
      break;
    case options::self:
      printf("# time   Fsk'   Fsk''\n");
      break;
    }
  }

  switch(opt.ave_type) {
  case options::numerical:
    for (int i=0; i<Fk.size(); i++)
      printf("%g %g %g\n",i*deltat,Fk[i].real(),Fk[i].imag());
  }
}

<<Fk.cc options>>=
class options {
public:
  enum   {numerical,analytical} ave_type;
  enum   {full,self} which_Fk;
  double k;
  int    nave;
  long   seed;
  bool   gromacs_trajectory;
  char **ifiles;

  options(int argc,char *argv[]);
  void show_usage(const char*);
} ;

options::options(int argc,char *argv[]) :
  ave_type(numerical),
  which_Fk(full),
  k(1),
  nave(1),
  seed(199332),
  gromacs_trajectory(false),
  ifiles(0)
{
  int c;
  while ((c=getopt(argc,argv,"FGShn:s:"))!=-1) {
    switch(c) {
    case 'F':
      which_Fk=full;
      break;
    case 'G':
      gromacs_trajectory=true;
      break;
    case 'S':
      which_Fk=self;
      break;
    case 'n':
      nave=atoi(optarg);
      break;
    case 's':
      seed=atol(optarg);
      break;
    case 'h':
    default:
      show_usage(argv[0]);
      break;
    }
  }
  if (argc-optind<2) show_usage(basename(argv[0]));
  k=atof(argv[optind]);
  ifiles=argv+optind+1;
} 

void options::show_usage(const char* name)
{
  std::cerr 
    << "usage: " << name << "[options] k file [file ...]\n\n"
    << "Computes the intermediate scattering functions F(k) or F_s(k) at the given\n"
    << "wavevector and from the given trajectory files.\n\n"
    << " Options:\n"
    << "   -h      show this help\n"
    << "   -F      compute F(k) [default]\n"
    << "   -S      compute only the self part F_s(k)\n"
    << "   -n nn   do nn averages over random directions of the wavevector\n"
    << "           (only applies to -F, with -S the average is over particles)\n"
    << "   -s xx   random number seed\n"
    << "   -G      read GROMACS .gro trajectories\n"
    << "\n";
  exit(1);
}



#endif 
