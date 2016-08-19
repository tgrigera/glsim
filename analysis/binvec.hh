/*
 * binvec.hh -- template class Binned_vector
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
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

#ifndef BINVEC_HH
#define BINVEC_HH

namespace glsim {

/** \class Binned_vector
    \ingroup Analysis

This is a vector that can be indexed with real (type double)
subindices.  It works like a histogram, with a range of subindices
mapping to a single bin, except that the bin can hold a value of any
type that is not necessarily an event count.

The number of bins and range is fixed on creation.  Bins are accessed
through operator[], which is overloaded to take integer values
(indexing bins) or double values (which will map to bins).  As in
std::vector, no bound checks are made.  If you need bound checks, use
at().

*/
template <typename T>
class Binned_vector {
public:
  Binned_vector(int nbins,double min,double max);

  /// Return number of bins
  int    nbins() const {return nbin;}
  double min() const {return min_;}
  double max() const {return max_;}
  /// Return bin width
  double delta() const {return delta_;}
  /// Return the center of the range corresponding to given bin
  double binc(int bin) {return min_+(bin+.5)*delta_;}
  /// Return the bin number corresponding to value r
  std::size_t binn(double r)
  {return (std::size_t) floor((r-min_)/delta_);}
  std::size_t binn_protected(double r);
  /// Element access without bound check
  T& operator[](std::size_t pos) {return data_[pos];}
  T& operator[](int pos) {return data_[(std::size_t) pos];}
  T& operator[](double r) {return data_[binn(r)];}
  /// Element access with bound check
  T& at(std::size_t pos);
  T& at(int pos) {return at((std::size_t) pos);}
  T& at(double r) {return data_[binn_protected(r)];}

private:
  long         out_below,out_above;
  double       min_,max_,delta_;
  long         nbin;
  long         ndata;
  std::valarray<T> data_;
} ;

template <typename T>
inline Binned_vector<T>::Binned_vector(int nbins,double min,double max) : 
  min_(min),
  max_(max),
  nbin(nbins),
  data_(nbins)
{
  delta_=(max_-min_)/nbins;
}

template <typename T>
std::size_t Binned_vector<T>::binn_protected(double r)
{
  std::size_t bin=floor((r-min_)/delta_);
  if (bin<0 || bin>=nbin) throw glsim::Out_of_range(HERE);
  return bin;
}

template <typename T>
T& Binned_vector<T>::at(std::size_t pos)
{
  if (pos<0 || pos>=nbin) throw glsim::Out_of_range(HERE);
  return data_[pos];
}

} /* namespace */


#endif /* BINVEC_HH */

