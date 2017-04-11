/*
 * gr.hh -- Classes for radial distribution function -- header
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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
