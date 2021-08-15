/*
 * mvector.hh -- Math vectors (header)
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019, 2020, 2021
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

#ifndef MVECTOR_HH
#define MVECTOR_HH

#include <math.h>
#include <iostream>
#include <array>

#include "exception.hh"
#include "cerrors.h"

namespace glsim {

/** \class mvector
 *  \ingroup Vector
 *  \brief Class for vector arithmetic
 *
 * This provides a class to represent a fixed-dimensional vector that
 * knows about vector operations like addition, dot product, etc.
 *
 * Implemented as a wrapper around std::array.  It is a wrapper and
 * not a simple template<> using so that operatior+= etc, can be
 * added.
 *
 */
template <unsigned int dimension>
class mvector {
public:
  std::array<double,dimension> data;

  /// \name Size
  /// @{
  size_t size() const noexcept
  {return data.size();}
  size_t dim() const noexcept
  {return data.size();}
  /// @}

  /// \name Access
  /// @{
  double& operator[](const size_t &i)
  {return data[i];}  ///< Acess element
  const double& operator[](const size_t &i) const
  {return data[i];}
  double& at(size_t i)
  {return data.at(i);}  ///< Acess element with bound check
  const double& at(size_t i) const
  {return data.at(i);}
  /// @}

  // Fill with value
  void fill(double &d)
  {data.fill(d);}   ///< Fill with value

  /// \name Addition and substraction
  /// @{
  mvector<dimension> &operator+=(const mvector<dimension>&);
  mvector<dimension> &operator-=(const mvector<dimension>&);
  const mvector<dimension> operator+(const mvector<dimension> &v) const
  { return mvector<dimension>(*this)+=v; }
  const mvector<dimension> operator-(const mvector<dimension> &v) const
  { return mvector<dimension>(*this)-=v; }
  /// @}

  /// \name Scalar product and multiplication by scalar
  /// @{
  mvector<dimension> &operator*=(const double &a);
  const mvector<dimension> operator*(const double &a) const
  { return mvector<dimension>(*this)*=a; }
  mvector<dimension> &operator/=(const double &a)
  { return operator*=(1./a); }
  double dot(const mvector<dimension>& v) const;
  /// @}

  /// \name Modulus, norm
  /// @{
  /// Squared modulus
  double modsq() const
  {return dot(*this);}
  double modulus() const
  {return sqrt(modsq());}
  /// Norm, same as modulus
  double norm() const
  {return modulus();}
  mvector<dimension>& normalize()
  { operator/=(norm()); return *this; }
  /// @}
} ;

typedef mvector<2> mvector_2D;

typedef mvector<3> mvector_3D;

// Implement addition and substraction

template <unsigned int dimension>
mvector<dimension> &mvector<dimension>::operator+=(const mvector<dimension>& v)
{
  for (int i=0; i<data.size(); ++i) data[i]+=v.data[i];
  return *this;
}

/// Template specializaiton for 2D addition
template<>
mvector<2>& mvector<2>::operator+=(const mvector<2>& v)
{
  data[0]+=v.data[0];
  data[1]+=v.data[1];
  return *this;
}

template<>
mvector<3>& mvector<3>::operator+=(const mvector<3>& v)
{
  data[0]+=v.data[0];
  data[1]+=v.data[1];
  data[2]+=v.data[2];
  return *this;
}

template <unsigned int dimension>
mvector<dimension> &mvector<dimension>::operator-=(const mvector<dimension>& v)
{
  for (int i=0; i<data.size(); ++i) data[i]-=v.data[i];
  return *this;
}

template<>
inline mvector<2>& mvector<2>::operator-=(const mvector<2>& v)
{
  data[0]-=v.data[0];
  data[1]-=v.data[1];
  return *this;
}

template<> inline
mvector<3>& mvector<3>::operator-=(const mvector<3>& v)
{
  data[0]-=v.data[0];
  data[1]-=v.data[1];
  data[2]-=v.data[2];
  return *this;
}

// Implement scalar product and multiplication by scalar

template <unsigned int dimension>
mvector<dimension> &mvector<dimension>::operator*=(const double &a)
{
  for (int i=0; i<data.size(); ++i) data[i]*=a;
  return *this;
}

template <> inline
mvector<2> &mvector<2>::operator*=(const double &a)
{
  data[0]*=a;
  data[1]*=a;
  return *this;
}

template <> inline
mvector<3> &mvector<3>::operator*=(const double &a)
{
  data[0]*=a;
  data[1]*=a;
  data[2]*=a;
  return *this;
}

template <unsigned int dimension>
const mvector<dimension> operator*(const double &a,const mvector<dimension>& v)
{
  return mvector<dimension>(v)*=a;
}

template <unsigned int dimension>
double mvector<dimension>::dot(const mvector<dimension>& v) const
{
  double d=data[0]*v.data[0];
  for (int i=1; i<data.size(); ++i) d+=data[i]*v.data[i];
  return d;
}

template <> inline
double mvector<2>::dot(const mvector<2>& v) const
{
  return data[0]*v.data[0] + data[1]*v.data[1];
}

template <> inline
double mvector<3>::dot(const mvector<3>& v) const
{
  return data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2];
}

// Vecctor product

const mvector_3D cross(const mvector_3D& a,const mvector_3D& b)
{
  return mvector_3D({a[1]*b[2] - a[2]*b[1],
		    a[2]*b[0] - a[0]*b[2],
		    a[0]*b[1] - a[1]*b[0]});
}


// I/O

/// \ingroup Vector
/// \brief Stream output for mvector
template <unsigned int dim>
std::ostream& operator<<(std::ostream& o,const mvector<dim> &v)
{
  o << '(';
  for (int i=0; i<v.dim()-1; ++i) o << v[i] << ", ";
  o << v[v.dim()-1] << ')';
  return o;
}

std::ostream& operator<<(std::ostream& o,const mvector_2D &v)
{
  o << '(' << v[0] << ", " << v[1] << ')';
  return o;
}

std::ostream& operator<<(std::ostream& o,const mvector_3D &v)
{
  o << '(' << v[0] << ", " << v[1] << ", " << v[2] << ')';
  return o;
}



} /* namespace */

#endif /* MVECTOR_HH */
