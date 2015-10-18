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

/** \class Fk
    \ingroup Structure
    \brief Intermediate scattering function

This class computes the intermediate scattering function, \f$F(k,t)\f$
for an isotropic system. Its definition is
\f[ F(k,t) = \frac{1}{N}\langle \rho_k(t)\rho_{-k}(0)\rangle =
\frac{1}{N} \sum_{i,j}\langle e^{i k\cdot[r_i(t)-r_j(0)]}. \f]

On construction the user sets the modulus of the wavevector, $k$, and
the class will choose a user-specified number of random directions and
average over them.  The calculatio uses the FFT-based time correlation
functions.

The way to use this class is to create an object, then push() the
desired configurations (assumed sequential and uniformly separated in
time) and finally call compute_Fk() to obtain the scattering function.
Data can be accessed through Fk_data() or printed in table form with
the overload of operator<<().
*/
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


/** \class Fsk
    \brief Self part of the intermediate scattering function

This class computes the self part of intermediate scattering function,
\f[ F_s(k,t) = \frac{1}{N} \sum_{i}\langle e^{i k\cdot[r_i(t)-r_i(0)]}
= \langle e^{-i k \cdot r_i(t)} e^{i k\cdot r_i(0)}.\f]

The modulus of the wavevector, $k$, is given on construction.  It is
computed along a random direction, and averaged over all particles.
The calculation uses the FFT-based time correlation functions.

The way to use this class is to create an object, then push() the
desired configurations (assumed sequential and uniformly separated in
time) and finally call compute_Fsk() to obtain the desired function.
Data can be accessed through Fsk_data() or printed in table form with
the overload of operator<<().
*/
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
