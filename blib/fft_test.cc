/*
 * fft_test.cc -- testing for FFT interface classes
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


#include <iostream>
#include <iomanip>
// #include <iterator>

#include "config.h"
#include "fft.hh"
#include "test_exception.hh"

void check_diff(const char* testing,double actual,double max)
{
  std::cout << "Testing " << testing << ": ";
  if (actual<max) {
    std::cout << "OK\n";
    return;
  }
  std::string exp=std::to_string(max);
  std::string act=std::to_string(actual);
  throw Test_failure(exp,act,HERE);
}


double msddiff(const glsim::vcomplex& a,const glsim::vcomplex& b)
{
  double d=0;
  for (int i=0; i<a.size(); ++i)
    d+=norm(a[i]-b[i]);
  d/=a.size();
  return d;
}

void show_test_result(int tst,const glsim::vcomplex& testd,
		      const glsim::vcomplex& testr,
		      const glsim::vcomplex& ft,const glsim::vcomplex& ift)
{
  double tol=1e-5;

  std::cout << "Test " << tst << "\n\n";
  std::cout << "Data \t\t\t\tFFT+Inverse\n";
  std::cout << std::scientific << std::setprecision(5) << std::showpos;
  for (int i=0; i<testd.size(); i++)
    std::cout << testd[i] << "\t" << ift[i] << '\n';
  std::cout << "\n" << std::noshowpos;
  double err=msddiff(testd,ift);
  std::cout << "ERR (MSD) = " << err << "\n\n";
  check_diff("FFT+Inverse",err,tol);

  std::cout << "ExpectedFFT\t\t\tActual-FFT\n" << std::showpos;
  for (int i=0; i<testd.size(); i++)
    std::cout << testr[i] << "\t" << ft[i] << "\n";
  std::cout << '\n' << std::noshowpos;
  err=msddiff(testr,ft);
  std::cout << "ERR (MSD) = " << err << "\n\n\n";
  check_diff("FFT+Expected",err,tol);
}

void copy_test_data(glsim::ComplexFFT &fft,glsim::vcomplex data)
{
  if (fft.tdata().size()!=data.size()) 
    fft.tdata_rw()=data;
  else for (int i=0; i<data.size(); i++) fft.tdatum(i)=data[i];
}

void load_test_data(int which,glsim::vcomplex& test,glsim::vcomplex& result)
{
  switch(which) {
  case 1:
    test.resize(8);
    test[0]=glsim::dcomplex(1,0);
    test[1]=glsim::dcomplex(0,0);
    std::copy(test.begin()+1,test.end()-1,test.begin()+2);
    result.resize(8);
    result[0]=glsim::dcomplex(8*0.125,0);
    std::copy(result.begin(),result.end()-1,result.begin()+1);
    break;
  case 2:
    test.resize(8);
    test[0]=glsim::dcomplex(0,0);
    test[1]=glsim::dcomplex(1,0);
    test[2]=glsim::dcomplex(0,0);
    std::copy(test.begin()+2,test.end()-1,test.begin()+3);
    // result should have constant modulus (equal to 1?)
    result.resize(8);
    result[0]=glsim::dcomplex(8,0)*glsim::dcomplex(+0.125,+0.000);
    result[1]=glsim::dcomplex(8,0)*glsim::dcomplex(+0.088,-0.088);
    result[2]=glsim::dcomplex(8,0)*glsim::dcomplex(+0.000,-0.125);
    result[3]=glsim::dcomplex(8,0)*glsim::dcomplex(-0.088,-0.088);
    result[4]=glsim::dcomplex(8,0)*glsim::dcomplex(-0.125,+0.000);
    result[5]=glsim::dcomplex(8,0)*glsim::dcomplex(-0.088,+0.088);
    result[6]=glsim::dcomplex(8,0)*glsim::dcomplex(+0.000,+0.125);
    result[7]=glsim::dcomplex(8,0)*glsim::dcomplex(+0.088,+0.088);
    break;
  }
}


void run_tests()
{
#ifdef HAVE_LIBFFTW3
  glsim::ComplexFFTW fft(glsim::FFT::out_of_place);
  std::cout << "FFT test\n\nLibrary: FFTW\n\n";
#else
  glsim::ComplexFFT_gsl_2n fft(glsim::FFT::out_of_place);
  std::cout << "FFT test\n\nLibrary: GSL\n\n";
#endif
  glsim::vcomplex    testd,testr;

  load_test_data(1,testd,testr);
  copy_test_data(fft,testd);
  fft.forward();
  fft.inverse();
  show_test_result(1,testd,testr,fft.fdata(),fft.tdata());
  load_test_data(2,testd,testr);
  copy_test_data(fft,testd);
  fft.forward();
  fft.inverse();
  show_test_result(2,testd,testr,fft.fdata(),fft.tdata());
}

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARD_ERROR 99

int main(int argc, char *argv[])
{
  int result=TEST_SUCCESS;
    
  try {
    run_tests();
  } catch(const Test_failure &e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
    result=TEST_FAILURE;
  } catch(const Test_skip &e) {
    std::cerr << e.what() << '\n';
    result=TEST_SKIPPED;
  } catch (const glsim::Runtime_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
    result=TEST_HARD_ERROR;
  } catch (const glsim::Logic_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
    result=TEST_HARD_ERROR;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << '\n';
    result=TEST_HARD_ERROR;
  }

  return result;
}
