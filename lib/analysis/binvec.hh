/*
 * binvec.hh -- template class Binned_vector
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

#ifndef BINVEC_HH
#define BINVEC_HH

#include <valarray>

#include "exception.hh"

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
  size_t size() const {return nbin;}
  double min() const {return min_;}
  double max() const {return max_;}
  /// Return bin width
  double delta() const {return delta_;}
  /// Return the center of the range corresponding to given bin
  double binc(int bin) const {return min_+(bin+.5)*delta_;}
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

  Binned_vector& operator=(const T& x)
  {data_=x; return *this;}

private:
  double       min_,max_,delta_;
  long         nbin;
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



/** \class Gird2D
    \ingroup Analysis

This is a two-dimensional array that can be indexed with real (type double)
subindices.  It is the 2-D version of Binned_vector.

The number of bins (or grid sites) and range is fixed on creation.
Bins are accessed through operator(), which is overloaded to take
integer values (selection by grid position) or double values (2-D
positions which will map to grid sites).  As in std::vector, no bound
checks are made.  If you need bound checks, use at().

*/
template <typename T>
class Grid2D {
public:
  Grid2D(double Lx,double Ly,double gsx,double gsy);

  /// Return number of bins
  int    nxbins() const {return nxbin;}
  int    nybins() const {return nybin;}
  size_t size() const {return nxbin*nybin;}
  double Lx() const {return Lx_;}
  double Ly() const {return Ly_;}
  /// Return grid size
  double gsx() const {return gsx_;}
  double gsy() const {return gsy_;}
  /// Return the center of the grid corresponding to given grid number
  void gridc(std::size_t gridn,double &xc,double &yc);
  /// Return the grid number corresponding to (x,y)
  std::size_t gridn(double x,double y) const
  {return (std::size_t) igrid(floor(x/gsx_),floor(y/gsy_));}
  std::size_t gridn_protected(double x,double y) const;
  /// Element access without bound check
  T& operator[](std::size_t gridn) {return data[gridn];}
  T& operator()(int nx,int ny) {return data[igrid(nx,ny)];}
  T& operator()(double px,double py) {return data[gridn(px,py)];}
  /// Element access with bound check
  T& at(int nx,int ny);
  T& at(double px,double py) {return data[gridn_protected(px,py)];}

private:
  double           Lx_,Ly_,gsx_,gsy_;
  double           nxbin,nybin;
  std::valarray<T> data;

  std::size_t igrid(int nx,int ny) const
  {return nxbin*ny + nx;}
  void gridnxny(int igrid,int &nx,int &ny)
  {ny = igrid/nxbin; nx=igrid - ny*nxbin;}
} ;

template <typename T>
inline Grid2D<T>::Grid2D(double Lx,double Ly,double gsx,double gsy) :
  Lx_(Lx), Ly_(Ly), gsx_(gsx), gsy_(gsy)
{
  nxbin=ceil(Lx/gsx_);
  nybin=ceil(Ly/gsy_);
  data.resize(nxbin*nybin);
}

template <typename T>
std::size_t Grid2D<T>::gridn_protected(double x,double y) const
{
  std::size_t igr=gridn(x,y);
  if (igr<0 || igr>=this->size()) throw glsim::Out_of_range(HERE);
  return igr;
}

template <typename T>
void Grid2D<T>::gridc(std::size_t gridn,double &xc,double &yc)
{
  int nx,ny;
  gridnxny(gridn,nx,ny);
  xc=(nx+0.5)*gsx_;
  yc=(ny+0.5)*gsy_;
}

template <typename T>
T& Grid2D<T>::at(int nx,int ny)
{
  if (nx<0 || nx>=nxbin) throw glsim::Out_of_range(HERE);
  if (ny<0 || ny>=nybin) throw glsim::Out_of_range(HERE);
  return data[igrid(nx,ny)];
}

  
} /* namespace */


#endif /* BINVEC_HH */

