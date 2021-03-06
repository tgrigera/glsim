/*
 * histogram.hh -- class for histogram generation
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

#ifndef HISTOGRAM_HH
#define HISTOGRAM_HH

#include <iostream>
#include <valarray>
#include <math.h>

namespace glsim {

class Histogram;

std::ostream &operator<<(std::ostream &os,Histogram &h);

/** \class Histogram
    \ingroup Analysis

This creates and holds an histogram of all the data added through the
push() method.  The histogram can be printed at any time through the
stream insertion operator overload, or individual bin counts examined.
The class keeps count of the outlier points (above and below the bin
range).

Output can be a count or a probability estimate.  Probabilities are
estimated taking outliers into account, so that the sum of the
probabilities will not necesarily equal 1, but rather will be \f$
N_\mathrm{tot}- N_\mathrm{outliers})/N_\mathrm{tot} \f$.
*/
class Histogram {
public:
  Histogram(int nbins,double min,double max,int extrabins=0);

  /// Set output mode (count or probability)
  Histogram& probability_output(bool n=true) {normalized_=n; return *this;}

  /// Add point
  int    push(double datum);
  /// Return number of bins
  int    nbins() const {return nbin;}
  /// Return number of points added so far
  long   npoints() const {return ndata;}
  double min() const {return min_;}
  double max() const {return max_;}
  /// Return bin width
  double delta() const {return delta_;}
  /// Return number of outliers
  long   outliers() const {return out_below+out_above;}
  /// Return number of points with value higher than max
  long   outliers_high() const {return out_above;}
  /// Return number of points with value lower than min
  long   outliers_low() const {return out_below;}
  /// Return the center of the range corresponding to given bin
  double binc(int bin) {return min_+(bin+.5)*delta_;}
  /// Return count value for the given bin
  long   count(int bin) const {return count_[bin];}
  /// Return the probability value for the given bin
  double prob(int bin) const {return (double) count_[bin]/(ndata*delta_);}
  /// Return histogram area (not counting outliers)
  double area();
  /// Compute the coarse-grained median (i.e. using bins and counts, not full list of data)
  double median();

private:
  long         out_below,out_above;
  double       min_,max_,delta_;
  long         nbin;
  long         ndata;
  std::valarray<long> count_;
  bool         normalized_;

  friend std::ostream &operator<<(std::ostream &os,Histogram &h);
} ;


/**
  Create histogram

 \param nbins    Number of bins
 \param min,max  Range
 \param extrabins Add extrabins bins below and above min and max, respectively, but
                  computing the bin width as (max-min)/nbins
*/
inline Histogram::Histogram(int nbins,double min,double max,int extrabins) : 
  min_(min),
  max_(max),
  nbin(nbins+2*extrabins),      // provide a padding above and below max/min
  ndata(0),
  count_(0L,nbins+2*extrabins),
  out_below(0), out_above(0),
  normalized_(false)
{
  delta_=(max_-min_)/nbins;
  max_+=extrabins*delta_;
  min_-=extrabins*delta_;
}

inline int Histogram::push(double datum) {
  ndata++;
  int bin= (int) floor((datum-min_)/delta_);
  if (bin<0) {out_below++; return -1;}
  else if (bin>=nbin) {out_above++; return -1;}
  else {count_[bin]++; return bin;}
}

inline double Histogram::area()
{
  return (double) (ndata-outliers())/ndata;
}

inline double Histogram::median()
{
  long cc=out_below;
  int  i;
  for (i=0; cc<npoints()/2; ++i)
    cc+=count(i);
  return binc(i);
}

/// Print histogram (two columns)
inline std::ostream &operator<<(std::ostream &os,Histogram &h)
{
  if (h.normalized_) {
    double fac=1./(h.ndata*h.delta_);
    for (int i=0; i<h.nbin; i++)
      os << h.min_+(i+.5)*h.delta_ << ' ' << fac*h.count_[i] << '\n';
  } else {
    for (int i=0; i<h.nbin; i++)
      os << h.min_+(i+.5)*h.delta_ << ' ' << h.count_[i] << '\n';
  }
  return os;
}


} /* namespace */

#endif /* HISTOGRAM_HH */

