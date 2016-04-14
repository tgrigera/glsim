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

#include <gsl/gsl_sf_trig.h>
#include "Sk.hh"

namespace glsim {

void Sk::basic_init(double box_length[])
{
  // choose deltak
  deltak_[0]=2*M_PI/box_length[0];
  deltak_[1]=2*M_PI/box_length[1];
  deltak_[2]=2*M_PI/box_length[2];

  sfact.clear();
  sfact.resize(Nk,0);
  nsamp=0;
}

Sk::Sk(double box_length[],int Nk_) :
  nsamp(0),
  Nk(Nk_),
  kdir(0)
{
  basic_init(box_length);
}

// this is isotropic
void Sk::push(OLconfiguration &conf)
{
  nsamp++;
  sfact[0]=conf.N;
  for (int ik=1; ik<Nk; ik++) {
    double k=ik*deltak_[0]/M_PI; // Divide by PI because of GSL's definition 
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


/**
  \param[in] kdir  Direction of k (0=x, 1=y, 2=z)

  We only allow the axes as directions, because in these directions
  the periodic boundary conditions are correctly enforced
  automatically by simply computing \f$\exp[-i k\cdot r]\f$.
*/  
void Sk::push(OLconfiguration &conf,int kdir_)
{
  nsamp++;
  kdir=kdir_;
  double kvec[3],deltakvec[3]={0.,0.,0.};
  std::complex<double> uimag(0,1);

  deltakvec[kdir]=deltak_[kdir];
  for (int ik=0; ik<Nk; ik++) {
    kvec[0]=ik*deltakvec[0];
    kvec[1]=ik*deltakvec[1];
    kvec[2]=ik*deltakvec[2];

    std::complex<double> S=0;

    for (int i=0; i<conf.N; i++) {
        double kr=kvec[0]*conf.r[i][0] + kvec[1]*conf.r[i][1] +
          kvec[2]*conf.r[i][2];
        S+=exp(-uimag*kr);
    }
    double Sr=(S.real()*S.real()+S.imag()*S.imag())/conf.N;

    // iterative average
    double Q=Sr-sfact[ik];
    sfact[ik]+=Q/nsamp;
  }
}

std::ostream& operator<<(std::ostream& o,Sk& S)
{
  for (int i=0; i<S.size(); ++i)
    o << S.k(i) << ' ' << S.S(i) << '\n';

  return o;
}

} /* namespace */
