/*
 * gr.cc -- Classes for radial distribution function
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

#include "gr.hh"

namespace glsim {

gr::gr(const OLconfiguration &c,double Dr_) :
  Dr(Dr_),
  nconf(0)
{
  double boxmin;
  boxmin=c.box_length[0];
  if (c.box_length[1]<boxmin) boxmin=c.box_length[1];
  if (c.box_length[2]<boxmin) boxmin=c.box_length[2];
  int Nbin=(int) ceil(0.5*boxmin/Dr);
  ggr.resize(Nbin,0.);
  init(c);
}

gr::gr(const OLconfiguration &c,int Nbin) :
  nconf(0)
{
  double boxmin=c.box_length[0];
  if (c.box_length[1]<boxmin) boxmin=c.box_length[1];
  if (c.box_length[2]<boxmin) boxmin=c.box_length[2];
  Dr=0.5*boxmin/Nbin;
  ggr.resize(Nbin,0.);
  init(c);
}

void gr::init(const glsim::OLconfiguration &c)
{
  double vol;
  vol=c.box_length[0]*c.box_length[1]*c.box_length[2];
  rhon=c.N*c.N/vol;
}

void gr::push(const glsim::OLconfiguration &c)
{
 int i,j;
 double temp;

 nconf++;
 for (int i=0; i<c.N-1; i++)
   for (int j=i+1; j<c.N; j++) {
     double r=sqrt(c.distancesq(i,j));
     int bin=(int) floor(r/Dr);
     if (bin<ggr.size())
       ggr[bin]+=2;
   }
}

double gr::rdf(int i) const
{
  double cst=16.*atan(1.)/3;
  double rmin=i*Dr;
  double rmax=rmin+Dr;
  double vol=cst*(rmax*rmax*rmax-rmin*rmin*rmin);
  return ggr[i]/(nconf*vol*rhon);
}


std::ostream& operator<<(std::ostream& o,const gr& rdf)
{
  o << "#  r   g(r)\n";
  for (int i=0; i<rdf.size(); i++)
    o << rdf.r(i) << " " << rdf.rdf(i) << '\n';
  return o;
}

} /* namespace */
