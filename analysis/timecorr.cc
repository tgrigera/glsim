/*
 * timecorr.cc -- Functions for discretized time correlations
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

/** \file timecorr.cc
    \ingroup TCorr
    \brief Functions for (discretized) time correlations 
*/

#include <algorithm>
#include <numeric>

#include "timecorr.hh"

namespace glsim {

/*****************************************************************************
 *
 * Self-correlations of scalar real quantities
 *
 */

//
// Direct method (O(N^2))
//

void correlation_1d_tti_direct(const std::vector<double> &a,
			       std::vector<double> &corr,int nt)
{
  int i;
  int n=a.size();
  if (nt>n) nt=n;
  corr.clear();
  corr.resize(nt,0);

  for (i=0; i<n; i++) {
    double ai=a[i];
    int nl=std::min(nt+i,n);
    for (int j=i; j<nl; j++)
      corr[j-i]+=ai*a[j];
  }

  for (i=0; i<nt; i++)
    corr[i]/=n-i;
}

void correlation_connected_1d_tti_direct(const std::vector<double> &a,
					 std::vector<double> &corr,int nt,
					 bool normalize)
{
  int i;
  int n=a.size();
  if (nt>n) nt=n;
  corr.clear();
  corr.resize(nt,0);

  /* Compute average */
  double ave=0;
  for (i=0; i<n; i++)
    ave+=a[i];
  ave/=n;

  /* Correlation loop */
  for (i=0; i<n; i++) {
    double ai=a[i]-ave;
    int nl=std::min(nt+i,n);
    for (int j=i; j<nl; j++)
      corr[j-i]+=ai*(a[j]-ave);
  }

  /* Normalize and return */
  double f=1;
  if (normalize) f=corr[0]/n;
  for (i=0; i<nt; i++)
    corr[i]/=f*(n-i);
}

void correlation_1d_tw_direct(const std::vector<double> &a,
			      std::vector<double> &corr,int tw)
{
  int i;
  int n=a.size();
  corr.clear();
  corr.resize(n-tw,0);

  /* Compute corr */
  for (i=0; i<n-tw; i++)
    corr[i]=a[tw]*a[tw+i];
}

void correlation_connected_1d_tw_direct(const std::vector<double> &a,
					std::vector<double> &corr,int tw,
					bool normalize)
{
  int i;
  int n=a.size();
  corr.clear();
  corr.resize(n-tw,0);

  /* Compute average */
  double ave=0;
  for (i=0; i<n; i++)
    ave+=a[i];
  ave/=n;

  for (i=0; i<n-tw; i++)
    corr[i]=(a[tw]-ave)*(a[tw+i]-ave);

  if (normalize)
    for (i=n-tw; i>=0; i--) corr[i]/=corr[0];
}

//
// With FFT  (O(N logN))
//

static long do_corr_1d_tti_fft_unnorm(glsim::RealFFT &corr)
{
  long N=corr.tdata().size();
  long N2=2*N;
  // Pad with 0s
  corr.tdata_rw().resize(N2,0.);
  // FFT
  corr.forward();
  // Compute squared modulus in half-complex mode
  corr.fdatum(0)*=corr.fdatum(0);
  for (int i=1; i<N; i++) {
    corr.fdatum(i)=corr.fdatum(i)*corr.fdatum(i) + 
                   corr.fdatum(N2-i)*corr.fdatum(N2-i);
    corr.fdatum(N2-i)=0;
  }
  corr.fdatum(N)*=corr.fdatum(N);
  // Invert FFT (without 1/N factor)
  corr.backward();
  // Return N so that caller can normalize (must divide by 2*N*(N-1)
  return N;
}

void correlation_1d_tti_fft(const std::vector<double> &a,RealFFT &corr,
			    int nt)
{
  corr.tdata_rw()=a;
  long N=do_corr_1d_tti_fft_unnorm(corr);
  corr.tdata_rw().resize(nt);  // nt tells how many points to keep, since correlation points near N have very large error

  // Normalize
  double N2=2*N;
  for (int i=0; i<nt; i++)
    corr.tdatum(i)/=N2*(N-i);
}

void correlation_connected_1d_tti_fft(const std::vector<double> &a,
				      glsim::RealFFT &corr,int nt,
				      bool normalize)
{
  int i;

  // Compute average and substract
  double ave=std::accumulate(a.begin(),a.end(),0.)/a.size();
  corr.tdata_rw().clear();
  corr.tdata_rw().resize(a.size());
  for (i=0; i<corr.tdata().size(); i++)
    corr.tdatum(i)=a[i]-ave;

  // Compute normal correlation and normalize
  long N=do_corr_1d_tti_fft_unnorm(corr);
  corr.tdata_rw().resize(nt);
  double f=2*N;
  corr.tdatum(0)/=f*N;
  if (normalize) {
    f*=corr.tdatum(0);
    corr.tdatum(0)=1;
  }
  for (int i=1; i<nt; i++)
    corr.tdatum(i)/=f*(N-i);
}


/*****************************************************************************
 *
 * Self-correlations of scalar complex quantities
 *
 */

//
// Direct method (O(N^2))
//

void correlation_1d_tti_direct(const vcomplex &a,vcomplex &corr,int nt)
{
  int i;
  int n=a.size();
  if (nt>n) nt=n;
  corr.clear();
  corr.resize(nt,0);

  for (i=0; i<n; i++) {
    dcomplex ai=a[i];
    int nl=std::min(nt+i,n);
    for (int j=i; j<nl; j++)
      corr[j-i]+=conj(ai)*a[j];
  }

  for (i=0; i<nt; i++)
    corr[i]/=n-i;
}

void correlation_connected_1d_tti_direct(const vcomplex &a,
					 vcomplex &corr,int nt,
					 bool normalize)
{
  int i;
  int n=a.size();
  if (nt>n) nt=n;
  corr.clear();
  corr.resize(nt,0);

  /* Compute average */
  dcomplex ave(0);
  for (i=0; i<n; i++)
    ave+=a[i];
  ave/=n;

  /* Correlation loop */
  for (i=0; i<n; i++) {
    dcomplex ai=a[i]-ave;
    int nl=std::min(nt+i,n);
    for (int j=i; j<nl; j++)
      corr[j-i]+=conj(ai)*(a[j]-ave);
  }

  /* Normalize and return */
  dcomplex f=1;
  if (normalize) f=corr[0]/dcomplex(n);
  for (i=0; i<nt; i++)
    corr[i]/=f*dcomplex(n-i);
}

//
// With FFT  (O(N logN))
//

static long do_corr_1d_tti_fft_unnorm(glsim::ComplexFFT &corr)
{
  long N=corr.tdata().size();
  long N2=2*N;
  // Pad with 0s
  corr.tdata_rw().resize(N2,0.);
  // FFT
  corr.forward();
  // Compute squared modulus in half-complex mode
  for (int i=0; i<N2; i++)
    corr.fdatum(i)*=std::conj(corr.fdatum(i));
  // Invert FFT (without 1/N factor)
  corr.backward();
  // Return N so that caller can normalize (must divide by 2*N*(N-i)
  return N;
}

void correlation_1d_tti_fft(const vcomplex &a,ComplexFFT &corr,int nt)
{
  corr.tdata_rw()=a;
  long N=do_corr_1d_tti_fft_unnorm(corr);
  // nt tells how many points to keep, points near N have very large error
  corr.tdata_rw().resize(nt);
  // Normalize
  double N2=2*N;
  for (int i=0; i<nt; i++)
    corr.tdatum(i)/=N2*(N-i);
}

void correlation_connected_1d_tti_fft(vcomplex &a,ComplexFFT &corr,int nt,
				      bool normalize)
{
  // Compute average and substract
  dcomplex ave=0;
  for (auto ai: a) ave+=ai;
  ave/=a.size();
  corr.tdata_rw().clear();
  corr.tdata_rw().resize(a.size());
  for (int i=0; i<corr.tdata().size(); i++)
    corr.tdatum(i)=a[i]-ave;

  // Compute correlation
  long N=do_corr_1d_tti_fft_unnorm(corr);
  corr.tdata_rw().resize(nt);
  // Normalize
  dcomplex f=2*N;
  corr.tdatum(0)/=f*dcomplex(N,0);
  if (normalize) {
    f*=corr.tdatum(0);
    corr.tdatum(0)=1;
  }
  for (int i=1; i<nt; i++)
    corr.tdatum(i)/=f*dcomplex(N-i,0);
}

} /* namespace */
