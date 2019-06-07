/*
 * iscatt.cc -- Classes to compute the intermediate scattering function
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

#include "iscatt.hh"

/*****************************************************************************
 *
 * F(k,t)
 *
 */

namespace glsim {

Fk::Fk(double k_,double deltat_,int Nav_) :
  k(k_),
  deltat(deltat_),
  Nav(Nav_),
  Npart(0),
  kr(Nav_),
  rhok_(Nav_)
{
  // fill random directions
  Spherical3d_distribution sr;
  for (int i=0; i<Nav; i++) {
    kr[i].resize(3);
    sr(&(kr[i][0]));
    for (int j=0; j<3; j++) kr[i][j]*=k;
  }
}

dcomplex Fk::rho_k(double r[][3],double k[])
{
  dcomplex rk(0);

  for (int i=0; i<Npart; i++) {
    double kr= k[0]*r[i][0] + k[1]*r[i][1] + k[2]*r[i][2];
    rk+=exp(dcomplex(0,-1)*kr);
  }
  return rk;
}

Fk& Fk::push_config(double r[][3],int N)
{
  Npart=N;

  for (int i=0; i<Nav; i++)
    rhok_[i].push_back(rho_k(r,&(kr[i][0])));

  return *this;
}

Fk& Fk::compute_Fk()
{
  // compute Fk for each direction and average
  dcomplex fac=1./dcomplex(Npart*Nav,0);
  int Fklen=rhok_[0].size()/2;
  Fk_.clear();
  Fk_.resize(Fklen,dcomplex(0,0));
  cFFT FF(FFT::in_place);
  for (int i=0; i<Nav; i++) {
    correlation_1d_tti_fft(rhok_[i],FF,Fklen);
    for (int j=0; j<FF.tdata_rw().size(); j++) Fk_[j]+=fac*FF.tdatum(j);
  }
  return *this;
}

std::ostream& operator<<(std::ostream& o,const Fk& Fk_)
{
  o << "# time   Fk'   Fk''\n";
  for (int i=0; i<Fk_.Fk_.size(); i++)
    o << i*Fk_.deltat << "  " << Fk_.Fk_[i].real()
      << "  " << Fk_.Fk_[i].imag() << '\n';
  return o;
}

/*****************************************************************************
 *
 * F_s(k,t)
 *
 */

Fsk::Fsk(double k_,double deltat_) :
  k(k_),
  deltat(deltat_),
  Npart(0)
{
  // choose a random dir
  Spherical3d_distribution sr;
  kr.resize(3);
  sr(&(kr[0]));
  for (int j=0; j<3; j++) kr[j]*=k;
}

Fsk& Fsk::push_config(double r[][3],int N)
{
  Npart=N;
  expkr.resize(N);
  
  for (int i=0; i<N; i++)
    expkr[i].push_back( exp( dcomplex(0,-1)*
			     (kr[0]*r[i][0] + kr[1]*r[i][1] + 
			      kr[2]*r[i][2]) ) );
  
  return *this;
}

Fsk& Fsk::compute_Fsk()
{
  dcomplex fac=dcomplex(1./Npart,0);
  int Fklen=expkr[0].size()/2;
  Fsk_.clear();
  Fsk_.resize(Fklen);
  cFFT FF(FFT::in_place);
  for (int i=0; i<Npart; i++) {
    correlation_1d_tti_fft(expkr[i],FF,Fklen);
    for (int j=0; j<FF.tdata_rw().size(); j++) Fsk_[j]+=fac*FF.tdatum(j);
  }
  return *this;
}

std::ostream& operator<<(std::ostream& o,const Fsk& Fsk_)
{
  o << "# time   F_s'(k)   F_s''(k)\n";
  for (int i=0; i<Fsk_.Fsk_.size(); i++)
    o << i*Fsk_.deltat << "  " << Fsk_.Fsk_[i].real()
      << "  " << Fsk_.Fsk_[i].imag() << '\n';
  return o;
}

} /* namespace */
