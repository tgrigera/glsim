/*
 * fft.hh -- Interface to FFT packages
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



#ifndef FFT_HH
#define FFT_HH

#include <vector>
#include <complex>
#include <gsl/gsl_errno.h>

#include "exception.hh"


namespace glsim {

/*****************************************************************************
 *
 * Interface classes
 *
 */

/** \class   FFT
    \ingroup FFT
    \brief   FFT interface definition
*/
class FFT {
public:
  enum storage_type {in_place,out_of_place} ;

  FFT(storage_type s) : storage(s) {}
  virtual ~FFT() {}
  virtual void forward()=0;
  virtual void backward()=0;
  virtual void inverse()=0;

protected:
  storage_type storage;
} ;

typedef std::complex<double>  dcomplex;
typedef std::vector<dcomplex> vcomplex;


/** \class   ComplexFFT
    \ingroup FFT
    \brief   Interface for FFT of complex data
*/
class ComplexFFT : public FFT {
public:
  ComplexFFT(storage_type s,int nreserve=0);
  ~ComplexFFT();

  const vcomplex &tdata() const {return *tdomain;}
  const vcomplex &fdata() const {return *fdomain;}
  dcomplex &tdatum(int n) const {return (*tdomain)[n];}
  dcomplex &fdatum(int n) const {return (*fdomain)[n];}
  virtual vcomplex &tdata_rw();
  virtual vcomplex &fdata_rw();

protected:
  vcomplex *tdomain,*fdomain;
} ;


typedef std::vector<double> vdouble;

/** \class   RealFFT
    \ingroup FFT
    \brief   Interface for FFT of real data
*/
class RealFFT : public FFT {
public:
  RealFFT(storage_type s,int nreserve=0);
  ~RealFFT();

  const vdouble &tdata() const {return *tdomain;}
  const vdouble &fdata() const {return *fdomain;}
  double &tdatum(int n) const {return (*tdomain)[n];}
  double &fdatum(int n) const {return (*fdomain)[n];}
  virtual vdouble &tdata_rw();
  virtual vdouble &fdata_rw();

protected:
  vdouble *tdomain,*fdomain;
} ;

/*****************************************************************************
 *
 * Classes for GSL routines
 *
 */

/** \class   ComplexFFT_gsl_2n
    \ingroup FFT
    \brief   Complex FFT through the Gnu Scientific Library

This uses the GSL complex routines for data sizes of powers of two.

*WARNING:* this code makes use of the C++11 standard guarantee that
allocation of complexes is contiguous (real, imag), and that vectors
are likewise contiguous (this is also guaranteed in C++03).
*/
class ComplexFFT_gsl_2n : public ComplexFFT {
public:
  ComplexFFT_gsl_2n(storage_type s,int nreserve=0) :  ComplexFFT(s,nreserve) {}
  void forward();
  void backward();
  void inverse();
} ;

/** \class   RealFFT_gsl_2n
    \ingroup FFT
    \brief   Real FFT through the Gnu Scientific Library

For FFT of real data. The transform is stored in GSLs half-complex
format also used by FFTW in some variants of the transforms.

It uses GSL real/half-complex routines for powers of two. These versions
manage the half-complex data in the packed form (see GSL reference,
15.6).

*/
class RealFFT_gsl_2n : public RealFFT {
public:
  RealFFT_gsl_2n(storage_type s,int nreserve=0) :  RealFFT(s,nreserve) {}
  void forward();
  void backward();
  void inverse();
} ;


/*****************************************************************************
 *
 * Classes for FFTW3 routines
 *
 */

#ifdef HAVE_LIBFFTW3
#include <fftw3.h>


/** \class   ComplexFFTW
    \ingroup FFT
    \brief   Complex FFT through the FFTW3 library
*/
class ComplexFFTW : public ComplexFFT {
public:
  ComplexFFTW(storage_type s,int n=0);
  ~ComplexFFTW();
  void forward();
  void backward();
  void inverse();

  vcomplex& tdata_rw();
  vcomplex& fdata_rw();

private:
  void       plan(bool time_rules);
  fftw_plan  forward_plan,backward_plan;
  bool       valid_plans;
} ;

inline ComplexFFTW::ComplexFFTW(storage_type s,int n) :
  ComplexFFT(s,n),
  forward_plan(0),
  backward_plan(0),
  valid_plans(false)
{}


/** \class   ComplexFFTW
    \ingroup FFT
    \brief   Complex FFT through the FFTW3 library
*/
class RealFFTW : public RealFFT {
public:
  RealFFTW(storage_type s,int n=0);
  ~RealFFTW();
  void forward();
  void backward();
  void inverse();

  vdouble& tdata_rw();
  vdouble& fdata_rw();

private:
  void       plan(bool time_rules);
  fftw_plan  r2hc,hc2r;
  bool       valid_plans;
} ;

inline RealFFTW::RealFFTW(storage_type s,int n) :
  RealFFT(s,n),
  r2hc(0),
  hc2r(0),
  valid_plans(false)
{}

#endif /* HAVE_LIBFFTW3 */

}

#endif /* FFT_HH */
