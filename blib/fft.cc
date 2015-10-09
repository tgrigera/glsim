/*
 * fft.cc -- Interface to FFT packages
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

#include <config.h>

#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>

#ifdef HAVE_LIBFFTW3
#include <algorithm>
#include <fftw3.h>
#endif

#include "fft.hh"

namespace glsim {

/*****************************************************************************
 *
 * Interface classes
 *
 */

ComplexFFT::ComplexFFT(storage_type s,int n) :
  FFT(s)
{
  tdomain=new vcomplex(n);
  if (storage==in_place)
    fdomain=tdomain;
  else 
    fdomain=new vcomplex(n);
}

ComplexFFT::~ComplexFFT()
{
  delete tdomain;
  if (storage==out_of_place) delete fdomain;
}

vcomplex& ComplexFFT::tdata_rw()
{
  return *tdomain;
}

vcomplex& ComplexFFT::fdata_rw()
{
  return *fdomain;
}

RealFFT::RealFFT(storage_type s,int n) :
  FFT(s)
{
  tdomain=new vdouble(n);
  if (storage==in_place)
    fdomain=tdomain;
  else 
    fdomain=new vdouble(n);
}

RealFFT::~RealFFT()
{
  delete tdomain;
  if (storage==out_of_place) delete fdomain;
}

vdouble& RealFFT::tdata_rw()
{
  return *tdomain;
}

vdouble& RealFFT::fdata_rw()
{
  return *fdomain;
}


/*****************************************************************************
 *
 * FFT with GSL
 *
 */


//
// Complex
//

void ComplexFFT_gsl_2n::forward()
{
  (*fdomain)=(*tdomain);
  gsl_complex_packed_array d=reinterpret_cast<double*>(&(*fdomain)[0]);
  int ret=gsl_fft_complex_radix2_forward(d,1,fdomain->size());
  if (ret) throw GSL_error(ret,HERE);
}

void ComplexFFT_gsl_2n::backward()
{
  (*tdomain)=(*fdomain);
  gsl_complex_packed_array d=reinterpret_cast<double*>(&(*tdomain)[0]);
  int ret=gsl_fft_complex_radix2_backward(d,1,tdomain->size());
  if (ret) throw GSL_error(ret,HERE);
}

void ComplexFFT_gsl_2n::inverse()
{
  *tdomain=*fdomain;
  gsl_complex_packed_array d=reinterpret_cast<double*>(&(*tdomain)[0]);
  int ret=gsl_fft_complex_radix2_inverse(d,1,tdomain->size());
  if (ret) throw GSL_error(ret,HERE);
}

//
// Real
//

void RealFFT_gsl_2n::forward()
{
  *fdomain=*tdomain;
  int ret=gsl_fft_real_radix2_transform(&(*fdomain)[0],1,fdomain->size());
  if (ret) throw GSL_error(ret,HERE);
}

void RealFFT_gsl_2n::inverse()
{
  *tdomain=*fdomain;
  int ret=gsl_fft_halfcomplex_radix2_inverse(&(*tdomain)[0],1,tdomain->size());
  if (ret) throw GSL_error(ret,HERE);
}

void RealFFT_gsl_2n::backward()
{
  *tdomain=*fdomain;
  int ret=gsl_fft_halfcomplex_radix2_backward(&(*tdomain)[0],1,tdomain->size());
  if (ret) throw GSL_error(ret,HERE);
}


/*****************************************************************************
 *
 * FFT with FFTW3
 *
 */

#ifdef HAVE_LIBFFTW3

//
// Complex
//

ComplexFFTW::~ComplexFFTW()
{
  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(backward_plan);
}

void ComplexFFTW::plan(bool time_rules)
{
  if (storage==out_of_place)
      if (tdomain->size()!=fdomain->size())
	if (time_rules) 
	  fdomain->resize(tdomain->size());
	else
	  tdomain->resize(fdomain->size());

  typedef double (*pp)[2];
  pp t=reinterpret_cast<pp>(&((*tdomain)[0]));
  pp f=reinterpret_cast<pp>(&((*fdomain)[0]));
  forward_plan=fftw_plan_dft_1d(tdomain->size(),t,f,
				FFTW_FORWARD,FFTW_ESTIMATE);
  backward_plan=fftw_plan_dft_1d(tdomain->size(),f,t,
				 FFTW_BACKWARD,FFTW_ESTIMATE);
  valid_plans=true;
}

void ComplexFFTW::forward()
{
  if (!valid_plans) plan(true);
  fftw_execute(forward_plan);
}

void ComplexFFTW::backward()
{
  if (!valid_plans) plan(false);
  fftw_execute(backward_plan);
}

void ComplexFFTW::inverse()
{
  backward();
  std::transform(tdomain->begin(),tdomain->end(),tdomain->begin(),
                 std::bind1st(std::multiplies<dcomplex>(),
			      dcomplex(1./tdomain->size())));  
}

vcomplex& ComplexFFTW::tdata_rw()
{
  valid_plans=false;
  return *tdomain;
}

vcomplex& ComplexFFTW::fdata_rw()
{
  valid_plans=false;
  return *fdomain;
}

//
// Real
//

RealFFTW::~RealFFTW()
{
  fftw_destroy_plan(r2hc);
  fftw_destroy_plan(hc2r);
}

void RealFFTW::plan(bool time_rules)
{
  if (storage==out_of_place)
      if (tdomain->size()!=fdomain->size())
	if (time_rules) 
	  fdomain->resize(tdomain->size());
	else
	  tdomain->resize(fdomain->size());

  r2hc=fftw_plan_r2r_1d(tdomain->size(),&(*tdomain)[0],&(*fdomain)[0],
			FFTW_R2HC,FFTW_ESTIMATE);
  // Always plan in-place because out-of-place destroys input anyway
  hc2r=fftw_plan_r2r_1d(tdomain->size(),&(*tdomain)[0],&(*tdomain)[0],
			FFTW_HC2R,FFTW_ESTIMATE);
  valid_plans=true;
}

void RealFFTW::forward()
{
  if (!valid_plans) plan(true);
  fftw_execute(r2hc);
}

void RealFFTW::backward()
{
  if (!valid_plans) plan(false);
  *tdomain=*fdomain;
  fftw_execute(hc2r);
}

void RealFFTW::inverse()
{
  backward();
  std::transform(tdomain->begin(),tdomain->end(),tdomain->begin(),
                 std::bind1st(std::multiplies<double>(),1./tdomain->size()));  
}

vdouble& RealFFTW::tdata_rw()
{
  valid_plans=false;
  return *tdomain;
}

vdouble& RealFFTW::fdata_rw()
{
  valid_plans=false;
  return *fdomain;
}

#endif /* HAVE_LIBFFTW3 */

}
