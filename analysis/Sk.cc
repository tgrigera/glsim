/*
 * Sk.cc -- Static structure factor class
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

#include <math.h>
#include <algorithm>

#include <gsl/gsl_sf_trig.h>
#include "Sk.hh"


namespace glsim {

Sk::Sk(double box_length[],int Nx,int Ny,int Nz) :
  nsamp(0),Sk0(0),
  Nkx(Nx), Nky(Ny), Nkz(Nz),
  sfact(0)
{
  if (Ny<0) Nkz=Nky=Nkx;

  deltak_[0]=2*M_PI/box_length[0];
  deltak_[1]=2*M_PI/box_length[1];
  deltak_[2]=2*M_PI/box_length[2];
  kmax=sqrt( Nkx*Nkx*deltak_[0]*deltak_[0] + Nky*Nky*deltak_[1]*deltak_[1] +
	     Nkz*Nkz*deltak_[2]*deltak_[2] );

  sfact = new glsim::Binned_vector<double>(std::max(Nkx, std::max(Nky,Nkz) ),0.,kmax);
  nsamp = new glsim::Binned_vector<double>(sfact->nbins(),0.,kmax);
  for (int i=0; i<sfact->nbins(); ++i)
    (*sfact)[i]=(*nsamp)[i]=0;
}

Sk::~Sk()
{
  delete sfact;
  delete nsamp;
}


/**
  \param[in] conf     Configuration

  Since we use only vectors belonging to the reciprocal lattice,
  periodic boundary conditions automatically enforced automatically by
  simply computing \f$\exp[-i k\cdot r]\f$.
*/  
void Sk::push(OLconfiguration &conf)
{
  double k,kx,ky,kz;
  std::complex<double> uimag(0,1);

  Sk0=conf.N;

  for (int ikx=0; ikx<Nkx; ++ikx) {
    kx=ikx*deltak_[0];
    for (int iky=0; iky<Nky; ++iky) {
      ky=iky*deltak_[1];
      for (int ikz=0; ikz<Nkz; ++ikz) {
	kz=ikz*deltak_[2];
	k=sqrt(kx*kx + ky*ky + kz*kz);

	if (k==0.) continue;  // S(k=0) = N, we know that and output separetely

	std::complex<double> rhok=0;
	for (int i=0; i<conf.N; i++) {
	  double kr=kx*conf.r[i][0] + ky*conf.r[i][1] + kz*conf.r[i][2];
	  rhok+=exp(-uimag*kr);
	}
	double Sk=(rhok.real()*rhok.real()+rhok.imag()*rhok.imag())/conf.N;

	// iterative average
	(*nsamp)[k]++;
	double Q=Sk-(*sfact)[k];
	(*sfact)[k]+=Q/(*nsamp)[k];
      }
    }
  }
}

std::ostream& operator<<(std::ostream& o,Sk& S)
{
  o << 0 << ' ' << S.Sk0 << '\n';
  for (int i=0; i<S.size(); ++i)
    o << S.k(i) << ' ' << S.S(i) << '\n';
  return o;
}

/*****************************************************************************
 *
 * Skiso
 *
 */
  
Skiso::Skiso(double box_length[],int Nk_) :
  nsamp(0),
  Nk(Nk_),
  sfact(0)
{
  double bl=std::min( std::min(box_length[0],box_length[1]) , box_length[2]);
  deltak_=2*M_PI/bl;
  kmax=Nk*deltak_;

  sfact = new double[Nk];
}

Skiso::~Skiso()
{
  delete[] sfact;
}

void Skiso::push(OLconfiguration &conf)
{
  nsamp++;
  sfact[0]=conf.N;
  for (int ik=1; ik<Nk; ik++) {
    double k=ik*deltak_/M_PI; // Divide by PI because of GSL's definition 
                                 // of the sinc function
    double S=0;
    for (int i=0; i<conf.N-1; i++)
      for (int j=i+1; j<conf.N; j++) {
        double kr=k*sqrt(conf.distancesq(i,j));
        S+=gsl_sf_sinc(kr); // S+=sin(kr)/kr;
    }
    S=1+2*S/conf.N;

    // iterative average
    double Q=S-sfact[ik];
    sfact[ik]+=Q/nsamp;
  }
}

std::ostream& operator<<(std::ostream& o,const Skiso& S)
{
  for (int i=0; i<S.size(); ++i)
    o << S.k(i) << ' ' << S.S(i) << '\n';
  return o;
}


} /* namespace */
