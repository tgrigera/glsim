/*
 * Sk.hh -- Header for structure factor class
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

#ifndef SK_HH
#define SK_HH

#include <iostream>
#include <vector>
#include <complex>
#include <math.h>

#include "olconfiguration.hh"

namespace glsim {

/** \class Sk
    \ingroup Structure
    \brief Compute the static structure factor

    This computes the static structure factor,
    \f[ S(k) = \frac{1}{N}\langle \rho_k \rho_{-k} \rangle =
    \frac{1}{N}\langle | \rho_k |^2 \rangle = \frac{1}{N} \sum_{ij}
    \left\langle e^{-i k\cdot r_i} e^{-k\cdot r_j} \right\rangle \f]


    The constructor takes the box size and number of desired k values
    as arguments.  Once constructed, the desired configurations are
    passed one by one to one of push() methdos, and the structure
    factor is computed and added to a running average.  The results
    can be read one at a time with the methods k() and S(), or can be
    sent as a whole to a stream with the operator<<().

    There are two versions of push().  The anisotropic version takes a
    vector as argument (apart from the configuration).  The structure
    factor is then computed along the direction of this vector.  The
    isotropic version computes the structure factor analytically
    averaged over all possible orientations, assuming an isotropic
    system, i.e.

    \f[ S_\text{iso}(k) = \sum_{ij} \frac{\sin(k |r_i-r_j|)}{k|r_i-r_j|} \f]
    
*/
class Sk {
public:
  ///  Constructor: give box size and desired number of k values
  Sk(double box_length[],int Nk);

  int    size() const;
  double deltak() const;
  double k(int) const;
  double S(int) const;

  /// Add the (isotropic) S(k) for the given configuration to the running average (O(N^2))
  void   push(OLconfiguration &conf);
  /// Add the value of S(k) in the given directon for the given configuration to the running average (O(N))
  void   push(OLconfiguration &conf,double kdir[]);

private:
  int    Nk;
  double deltak_;

  std::vector<double> sfact;
  int                 nsamp;

  void basic_init(double box_length[]);

  friend std::ostream& operator<<(std::ostream&,Sk&);
} ;

/// Number of k (scattering vector) values
inline int    Sk::size() const {return Nk;}

/// Returs \f$\Delta k\f$
inline double Sk::deltak() const {return deltak_;}

/// Returns the (modulus of) the ith scattering vector
inline double Sk::k(int i) const {return i*deltak_;}

/// Returns the structure factor at the ith scattering vector
inline double Sk::S(int i) const {return sfact.at(i);}

/// Print the structure factor for all vectors in two columns
std::ostream& operator<<(std::ostream&,Sk&);

} /* namespace */

#endif /* SK_HH */

