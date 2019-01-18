/*
 * gr.hh -- Classes for radial distribution function -- header
 *
 * This file is part of olglsim, a numerical simulation class library and
 * helper programs.
 *
 * olglsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * olglsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use olglsim to produced published work, or if you redistribute a
 * modified version of glsim, or code based on glsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * olglsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding. 
 *
 * olglsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

#ifndef GR_HH
#define GR_HH

#include "olconfiguration.hh"

namespace glsim {

/** \class gr
    \ingroup Structure

  Compute the radial distribution function (\f$g(r)\f$)
*/
class gr {
public:
  /// Constructor: give desired \f$\Delta r\f$
  gr(const glsim::OLconfiguration& c,double Dr);
  /// Constructor: give desired number of r values
  gr(const glsim::OLconfiguration& c,int Nbin);

  /// Number of bins
  int    size() const {return ggr.size();}
  /// \f$\Delta r\f$
  double deltar() const {return Dr;}
  /// Radius corresponding to ith bin
  double r(int i) const {return (i+0.5)*Dr;}
  /// g(r) at ith bin
  double rdf(int i) const;

  void push(const OLconfiguration&);

private:
  void   init(const OLconfiguration&);

  int    nconf;
  double Dr;
  std::vector<double> ggr;
  double  rhon;
} ;

std::ostream& operator<<(std::ostream& o,const gr& rdf);

}

#endif /* GR_HH */
