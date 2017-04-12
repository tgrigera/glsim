/*
 * Grk.cc -- Space-correlation classes
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017
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
#include "Grk.hh"

namespace glsim {

Gk::Gk(double box_length[],int Nx,int Ny,int Nz) :
  nsamp(0), Gk0(0.), nsamp0(0),
  Nkx(Nx), Nky(Ny), Nkz(Nz),
  gk(0)
{
  if (Ny<0) Nkz=Nky=Nkx;

  deltak_[0]=2*M_PI/box_length[0];
  deltak_[1]=2*M_PI/box_length[1];
  deltak_[2]=2*M_PI/box_length[2];
  kmax=sqrt( Nkx*Nkx*deltak_[0]*deltak_[0] + Nky*Nky*deltak_[1]*deltak_[1] +
	     Nkz*Nkz*deltak_[2]*deltak_[2] );

  int nbins = std::max(Nkx, std::max(Nky,Nkz) );
  gk = new glsim::Binned_vector<double>(nbins,0.,kmax);
  nsamp = new glsim::Binned_vector<double>(gk->nbins(),0.,kmax);
  for (int i=0; i<gk->nbins(); ++i)
    (*gk)[i]=(*nsamp)[i]=0;
}

Gk::~Gk()
{
  delete gk;
  delete nsamp;
}

/**
  \param[in] conf     Configuration

  Since we use only vectors belonging to the reciprocal lattice,
  periodic boundary conditions are automatically enforced by
  simply computing \f$\exp[-i k\cdot r]\f$.
*/  
void Gk::push(OLconfiguration &conf,double phi[])
{
  double k,kx,ky,kz;
  std::complex<double> uimag(0,1);

  for (int ikx=0; ikx<Nkx; ++ikx) {
    kx=ikx*deltak_[0];
    for (int iky=0; iky<Nky; ++iky) {
      ky=iky*deltak_[1];
      for (int ikz=0; ikz<Nkz; ++ikz) {
	kz=ikz*deltak_[2];
	std::complex<double> rhok=0;
	for (int i=0; i<conf.N; i++) {
	  double kr=kx*conf.r[i][0] + ky*conf.r[i][1] + kz*conf.r[i][2];
	  rhok+=phi[i]*exp(-uimag*kr);
	}
	double Gk=norm(rhok)/conf.N;

	k=sqrt(kx*kx + ky*ky + kz*kz);
	// iterative average
	if (k==0.) {
	  nsamp0++;
	  double Q=Gk-Gk0;
	  Gk0+=Q/nsamp0;
	} else {
	  (*nsamp)[k]++;
	  double Q=Gk-(*gk)[k];
	  (*gk)[k]+=Q/(*nsamp)[k];
	}
      }
    }
  }
}

void Gk::push(OLconfiguration &conf,double phi[][3])
{
  double k,kx,ky,kz;
  std::complex<double> uimag(0,1);

  for (int ikx=0; ikx<Nkx; ++ikx) {
    kx=ikx*deltak_[0];
    for (int iky=0; iky<Nky; ++iky) {
      ky=iky*deltak_[1];
      for (int ikz=0; ikz<Nkz; ++ikz) {
	kz=ikz*deltak_[2];

	std::complex<double> rhokx=0,rhoky=0,rhokz=0;
	for (int i=0; i<conf.N; i++) {
	  double kr=kx*conf.r[i][0] + ky*conf.r[i][1] + kz*conf.r[i][2];
	  rhokx+= phi[i][0]*exp(-uimag*kr);
	  rhoky+= phi[i][1]*exp(-uimag*kr);
	  rhokz+= phi[i][2]*exp(-uimag*kr);
	}
	double Gk=(norm(rhokx) + norm(rhoky) + norm(rhokz) )/conf.N;

	k=sqrt(kx*kx + ky*ky + kz*kz);
	// iterative average
	if (k==0.) {
	  nsamp0++;
	  double Q=Gk-Gk0;
	  Gk0+=Q/nsamp0;
	} else {
	  (*nsamp)[k]++;
	  double Q=Gk-(*gk)[k];
	  (*gk)[k]+=Q/(*nsamp)[k];
	}
      }
    }
  }
}

std::ostream& operator<<(std::ostream& o,const Gk& G)
{
  o << 0 << ' ' << G.Gk0 << '\n';
  for (int i=0; i<G.size()-1; ++i) {
    if ((*G.nsamp)[i]>0)
      o << G.k(i) << ' ' << G.G(i) << '\n';
  }
  return o;
}


/******************************************************************************
 *
 * Gkcub
 *
 */
// Gkcub::Gkcub(double box_length[],int Nk_) :
//   nsamp(0), Nk(Nk_),
//   sfact(0)
// {
//   deltak_=2*M_PI/box_length[0];
//   nbins = 3*Nk*Nk;
//   sfact = new double[nbins];
//   nsamp = new long[nbins];
//   for (int i=0; i<nbins; ++i)
//     sfact[i]=nsamp[i]=0;
// }

// Gkcub::~Gkcub()
// {
//   delete[] sfact;
//   delete[] nsamp;
// }


// void Gkcub::push(OLconfiguration &conf)
// {
//   double k,kx,ky,kz;
//   std::complex<double> uimag(0,1);

//   sfact[0]=conf.N; nsamp[0]++;

//   for (int ikx=0; ikx<Nk; ++ikx) {
//     kx=ikx*deltak_;
//     for (int iky=0; iky<Nk; ++iky) {
//       ky=iky*deltak_;
//       for (int ikz=0; ikz<Nk; ++ikz) {
// 	kz=ikz*deltak_;
//         long ik=ikx*ikx+iky*iky+ikz*ikz;
// 	if (ik==0) continue;  // S(k=0) = N, we know that and output separetely

// 	std::complex<double> rhok=0;
// 	for (int i=0; i<conf.N; i++) {
// 	  double kr=kx*conf.r[i][0] + ky*conf.r[i][1] + kz*conf.r[i][2];
// 	  rhok+=exp(-uimag*kr);
// 	}
// 	double Gkcub=(rhok.real()*rhok.real()+rhok.imag()*rhok.imag())/conf.N;

// 	// iterative average
// 	nsamp[ik]++;
// 	double Q=Gkcub-sfact[ik];
// 	sfact[ik]+=Q/nsamp[ik];
//       }
//     }
//   }
// }

// std::ostream& operator<<(std::ostream& o,const Gkcub& S)
// {
//   for (int i=0; i<S.size()-1; ++i) {
//     if (S.nsamp[i]>0)
//       o << S.k(i) << ' ' << S.S(i) << '\n';
//   }
//   return o;
// }

/*****************************************************************************
 *
 * Grk
 *
 */
Grk::Grk(double box_length[],int Nr,int Nk_,bool space) :
  Gr_space(space), Nk(Nk_), gr(0), gk(0), nsr(0), nsk(0),
  gr0(0), nsr0(0)
{
  double rmax=0.51*sqrt(box_length[0]*box_length[0] + box_length[1]*box_length[1] +
		       box_length[2]*box_length[2]);
  vol=box_length[0]*box_length[1]*box_length[2];
  double bl=std::min( std::min(box_length[0],box_length[1]) , box_length[2]);
  deltak_=2*M_PI/bl;

  gr = new Binned_vector<double>(Nr,0.,rmax);
  for (int i=0; i<gr->size(); ++i) (*gr)[i]=0;
  if (!Gr_space) {
    nsr = new Binned_vector<long>(Nr,0.,rmax);
    for (int i=0; i<nsr->size(); ++i) (*nsr)[i]=0;
  }
  gk = new double[Nk];
  for (int i=0; i<Nk; ++i) gk[i]=0;
}

Grk::~Grk()
{
  delete gr,nsr;
  delete[] gk;
}

inline void Grk::iteav(double &ave,long &count,double &datum)
{
  ++count;
  double Q=datum-ave;
  ave+=Q/count;
}

void Grk::push(OLconfiguration &conf,double phi[])
{
  double r,k,sGr,sGk;
  
  N=conf.N;
  ++nsk;

  if (Gr_space) {
    for (int i=0; i<conf.N; ++i) {
      sGr=phi[i]*phi[i];
      gr0+=sGr;
      gk[0]+=sGr; // Because sinc=1 in this case (r=0)
      for (int ik=1; ik<Nk; ik++) {
	k=ik*deltak_/M_PI; // Divide by PI because of GSL's definition 
	gk[ik]+=sGr;
      }
    }
    for (int i=0; i<conf.N-1; i++) {
      for (int j=i+1; j<conf.N; j++) {
	r=sqrt(conf.distancesq(i,j));
	sGr=phi[i]*phi[j];
	(*gr)[r]+=2*sGr;
	gk[0]+=2*sGr; // Becuase sinc=1 in this case (k=0), and counts twice because i!=j
	for (int ik=1; ik<Nk; ik++) {
	  double k=ik*deltak_/M_PI; // Divide by PI because of GSL's definition
	  sGk=phi[i]*phi[j]*gsl_sf_sinc(k*r); // S+=phi_*phi_j*sin(kr)/kr;
	  gk[ik]+=2*sGk;
	}
      }
    }
  } else {  // Is the same as before but for G(r) the number of pairs is counted
    for (int i=0; i<conf.N; ++i) {
      sGr=phi[i]*phi[i];
      iteav(gr0,nsr0,sGr);
      gk[0]+=sGr; // Because sinc=1 in this case (r=0)
      for (int ik=1; ik<Nk; ik++) {
	k=ik*deltak_/M_PI; // Divide by PI because of GSL's definition 
	gk[ik]+=sGr;
      }
    }
    for (int i=0; i<conf.N-1; i++) {
      for (int j=i+1; j<conf.N; j++) {
	r=sqrt(conf.distancesq(i,j));
	sGr=phi[i]*phi[j];
	iteav((*gr)[r],(*nsr)[r],sGr);
	iteav((*gr)[r],(*nsr)[r],sGr);  // The pair is counted twice

	gk[0]+=2*sGr; // Becuase sinc=1 in this case (k=0), and counts twice because i!=j
	for (int ik=1; ik<Nk; ik++) {
	  double k=ik*deltak_/M_PI; // Divide by PI because of GSL's definition
	  sGk=phi[i]*phi[j]*gsl_sf_sinc(k*r); // S+=phi_*phi_j*sin(kr)/kr;
	  gk[ik]+=2*sGk;
	}
      }
    }
  }
}

void Grk::push(OLconfiguration &conf,double phi[][3])
{
  double r,k,sGr,sGk;
  
  N=conf.N;
  ++nsk;

  if (Gr_space) {
    for (int i=0; i<conf.N; ++i) {
      sGr=dotp(phi[i],phi[i]);
      gr0+=sGr;
      gk[0]+=sGr; // Because sinc=1 in this case (r=0)
      for (int ik=1; ik<Nk; ik++) {
	k=ik*deltak_/M_PI; // Divide by PI because of GSL's definition 
	gk[ik]+=sGr;
      }
    }
    for (int i=0; i<conf.N-1; i++) {
      for (int j=i+1; j<conf.N; j++) {
	r=sqrt(conf.distancesq(i,j));
	sGr=dotp(phi[i],phi[j]);
	(*gr)[r]+=2*sGr;
	gk[0]+=2*sGr; // Becuase sinc=1 in this case (k=0), and counts twice because i!=j
	for (int ik=1; ik<Nk; ik++) {
	  double k=ik*deltak_/M_PI; // Divide by PI because of GSL's definition
	  sGk=dotp(phi[i],phi[j])*gsl_sf_sinc(k*r); // S+=phi_*phi_j*sin(kr)/kr;
	  gk[ik]+=2*sGk;
	}
      }
    }
  } else {  // Is the same as before but for G(r) the number of pairs is counted
    for (int i=0; i<conf.N; ++i) {
      sGr=dotp(phi[i],phi[i]);
      iteav(gr0,nsr0,sGr);
      gk[0]+=sGr; // Because sinc=1 in this case (r=0)
      for (int ik=1; ik<Nk; ik++) {
	k=ik*deltak_/M_PI; // Divide by PI because of GSL's definition 
	gk[ik]+=sGr;
      }
    }
    for (int i=0; i<conf.N-1; i++) {
      for (int j=i+1; j<conf.N; j++) {
	r=sqrt(conf.distancesq(i,j));
	sGr=dotp(phi[i],phi[j]);
	iteav((*gr)[r],(*nsr)[r],sGr);
	iteav((*gr)[r],(*nsr)[r],sGr);  // The pair is counted twice

	gk[0]+=2*sGr; // Becuase sinc=1 in this case (k=0), and counts twice because i!=j
	for (int ik=1; ik<Nk; ik++) {
	  double k=ik*deltak_/M_PI; // Divide by PI because of GSL's definition
	  sGk=dotp(phi[i],phi[j])*gsl_sf_sinc(k*r); // S+=phi_*phi_j*sin(kr)/kr;
	  gk[ik]+=2*sGk;
	}
      }
    }
  }
}

double Grk::Gr(int i) const
{
  if (Gr_space) {
    double rhon=N*N/vol;
    double cst=4*M_PI/3;
    double rc=gr->binc(i);
    double rmin=rc-0.5*gr->delta();
    double rmax=rc+0.5*gr->delta();
    double vol=cst*(rmax*rmax*rmax-rmin*rmin*rmin);
    return (*gr)[i]/(nsk*vol*rhon);
  } else
    return (*gr)[i]/gr->delta();
}

std::ostream& operator<<(std::ostream& o,const Grk::prtGr Gp)
{
  const Grk *G=Gp.p;
  o << 0 << ' ' << G->gr0 << '\n';
  for (int i=0; i<G->sizer(); ++i)
    o << G->r(i) << ' ' << G->Gr(i) << '\n';
  return o;
}

std::ostream& operator<<(std::ostream& o,const Grk::prtGk Gp)
{
  const Grk *G=Gp.p;
  for (int i=0; i<G->sizek(); ++i)
    o << G->k(i) << ' ' << G->Gk(i) << '\n';
  return o;
}

} /* namespace */
