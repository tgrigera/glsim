/*
 * Grk.hh -- Header for space-correlatoins classes
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

#ifndef GRK_HH
#define GRK_HH

#include <iostream>
#include <vector>
#include <complex>
#include <math.h>

#include "olconfiguration.hh"
#include "binvec.hh"

namespace glsim {

/** \class Gk
    \ingroup Structure
    \brief Compute the Fourier space correlation function

    This computes the correlation function in \f$k\f$-space
    \f[
        G(k) = \frac{1}{N}\langle \rho_k \rho_{-k} \rangle =
       \frac{1}{N}\langle | \phi_k |^2 \rangle = \frac{1}{N} \sum_{ij}
       \left\langle \phi_i \phi_j  e^{-i k\cdot r_i} e^{-k\cdot r_j} \right\rangle 
    \f]

    This is like the structure factor computed by Sk but for property
    \f$\phi\f$ associated to each particle.

    The result is the correlation as a function of the modulus of k,
    obtained by averaging over all the compatible wavevectors that
    belong to the reciprocal lattice of the Bravais lattice
    corresponding to the periodic repetitions of the simulation box.

    The constructor takes the box size and number of desired k values
    as arguments.  Once constructed, the desired configurations are
    passed one by one the push() methd, and the structure factor is
    computed and added to a running average.  The results can be read
    one at a time with the methods k() and S(), or can be sent as a
    whole to a stream with the operator<<().

    The push() function is linear in the number of particle but cubic
    in Nx (or linear in Nx*Ny*Nz).
*/
class Gk {
public:
  ///  Constructor: give box size and desired number of k values on each axis (if only one given, the same number is used for all three)
  Gk(double box_length[],int Nx,int Ny=-1,int Nz=-1);
  ~Gk();

  int    size() const;
  double deltak() const;
  double k(int) const;
  double G(int) const;

  /// Add the value of S(k) in the given directon for the given configuration to the running average (O(N))
  void   push(OLconfiguration &conf,double phi[]);
  void   push(OLconfiguration &conf,double phi[][3]);

private:
  int    Nkx,Nky,Nkz;
  double deltak_[3];
  double kmax;

  double Gk0;
  double nsamp0;
  Binned_vector<double> *gk,*nsamp;

  friend std::ostream& operator<<(std::ostream& o,const Gk& S);
} ;

/// Number of k (scattering vector) values
inline int    Gk::size() const {return gk->nbins();}

/// Returns \f$\Delta k\f$.  This is not the value used for the vector
/// wavelength (which can be different for the different axes) but the
/// bin witdh used to average directions
inline double Gk::deltak() const {return gk->delta();}

/// Returns the (modulus of) the ith scattering vector
inline double Gk::k(int i) const {return gk->binc(i);}

/// Returns the correlation at the ith scattering vector
inline double Gk::G(int i) const {return gk->at(i);}

/// Print the correlation function for all wavevectors in two columns
std::ostream& operator<<(std::ostream&,const Gk&);


// \class Gkcub
//     \ingroup Structure
//     \brief Compute the static structure factor (for square boxes)

//     This is like Gk, but performs no binning and returns the actual
//     values of the wave vector used, without averaging close values of
//     k as Gk does.  Only good for cubic simulation boxes.
// class Gkcub {
// public:
//   ///  Constructor: give box size and desired number of k values on each axis (if only one given, the same number is used for all three)
//   Gkcub(double box_length[],int Nk);
//   ~Gkcub();

//   int    size() const;
//   double deltak() const;
//   double k(int) const;
//   double S(int) const;

//   /// Add the value of S(k) in the given directon for the given configuration to the running average (O(N))
//   void   push(OLconfiguration &conf);

// private:
//   int    Nk;
//   double deltak_;
//   double kmax;

//   double *sfact;
//   long   nbins,*nsamp;

//   friend std::ostream& operator<<(std::ostream& o,const Gkcub& S);
// } ;

// /// Number of k (scattering vector) values
// inline int  Gkcub::size() const {return nbins;}

// /// Returns \f$\Delta k\f$.  This is not the value used for the vector
// /// wavelength (which can be different for the different axes) but the
// /// bin witdh used to average directions
// inline double Gkcub::deltak() const {return deltak_;}

// /// Returns the (modulus of) the ith scattering vector
// inline double Gkcub::k(int i) const {return deltak_*sqrt(i); }

// /// Returns the structure factor at the ith scattering vector
// inline double Gkcub::S(int i) const {return sfact[i];}

// /// Print the structure factor for all vectors in two columns
// std::ostream& operator<<(std::ostream&,const Gkcub&);


/** \class Grk
    \ingroup Structure
    \brief Compute the space correlation function in real and Fourier space (isotropic)

    This computes the space correlation function of some property \f$\phi\f$,
    \f[ G(r) = \frac{ \sum_{ij} \phi_i\phi_j \delta (r-r_{ij})}{\sum_{ij} \delta(r-r_{ij})}, \f]
    and its (analytically averaged over directions) Fourier transform
    \f[ G(k) = \frac{1}{N} \sum_{ij} \phi_i \phi_j \frac{\sin kr_{ij}}{kr_{ij}}. \f]

    \f$\phi\f$ can be a scalar or a 3-\f$d\f$ vector, in this case the
    product is understood as a scalar product.

    Optionally, the real space correlation can be computed with the
    following alternative definition, based on space-averaged (rather
    than particle-averaged) fields:
    \f[ G^{(S)}(\mathbf{r})  =  \frac{1}{N \rho_0} \sum_{ij} \phi_i \phi_j \delta ({\bf r} -{\bf r}_{ij}), \f]
    where \f$\rho_0=N/V\f$ is the number density.

    Note that the sum includes the terms with
    \f$i=j\f$, so that \f$\lim_{k\to\infty}G(k)=\sum_i \phi^2_i /N\f$.

    This class is used in the same way as glsim::Gk.

    The push() function is linear in the number of wavevectors but
    quadratic in the number of particles.
*/
class Grk {
public:
  class prtGk;
  class prtGr; 

  ///  Constructor: give box size and desired number of k values
  Grk(double box_length[],int Nr,int Nk,bool Gr_space);
  ~Grk();

  int    sizer() const,sizek() const;
  double deltar() const,deltak() const;
  double r(int) const,k(int) const;
  /// Returns G(r) at the ith distance
  double Gr(int) const;
  double Gk(int) const;
  
  prtGk  pGk() const;
  prtGr  pGr() const;

  /// Update the running average of \f$G(r)\f$ and \f$G(k)\f$ using the given configuration and property phi
  void   push(OLconfiguration &conf,double phi[]);
  /// Update the running average of \f$G(r)\f$ and \f$G(k)\f$ using the given configuration and 3-d vector property  phi
  void   push(OLconfiguration &conf,double phi[][3]);

private:
  bool   Gr_space;
  void   iteav(double&,long&,double&);
  double dotp(double a[3],double b[3]);

  int    N,Nk;
  double vol,deltak_;
  Binned_vector<double> *gr;
  Binned_vector<long>   *nsr;
  double                gr0,*gk;
  long                  nsk,nsr0;

  friend std::ostream& operator<<(std::ostream&,const prtGk);
  friend std::ostream& operator<<(std::ostream&,const prtGr);

} ;

class Grk::prtGr {
public:
  prtGr(const Grk* g) : p(g) {}
  const Grk *p;
} ;

class Grk::prtGk {
public:
  prtGk(const Grk* g) : p(g) {}
  const Grk *p;
} ;

inline Grk::prtGr Grk::pGr() const {return prtGr(this);}

inline Grk::prtGk Grk::pGk() const {return prtGk(this);}

inline double Grk::dotp(double a[3],double b[3])
{  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

/// Number of r (space) values
inline int  Grk::sizer() const {return gr->nbins();}

/// Number of k (scattering vector) values
inline int  Grk::sizek() const {return Nk;}

/// Returns \f$\Delta r\f$
inline double Grk::deltar() const {return gr->delta();}

/// Returns \f$\Delta k\f$
inline double Grk::deltak() const {return deltak_;}

/// Returns the ith distance
inline double Grk::r(int i) const {return gr->binc(i);}

/// Returns the (modulus of) the ith scattering vector
inline double Grk::k(int i) const {return i*deltak_;}

/// Returns G(k) at the ith scattering vector
inline double Grk::Gk(int i) const {return gk[i]/(N*nsk);}

/// Print G(k) for all vectors in two columns
std::ostream& operator<<(std::ostream& o,const Grk::prtGr G);

/// Print G(k) for all vectors in two columns
std::ostream& operator<<(std::ostream& o,const Grk::prtGk G);


} /* namespace */

#endif /* GRK_HH */
