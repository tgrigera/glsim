/*
 * geoave.cc -- Class to average in geometrically growing windows
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

#include "geoave.hh"

namespace glsim {

Geoave::Geoave(double t0_,double wfactor_,double base_) :
  base(base_),
  t0(t0_),
  wfactor(wfactor_)
{
  logwf=log(wfactor);
  read_fb=(wfactor-1)/base;
  count.reserve(20);
  rave.reserve(20);
  rvarn.reserve(20);
}

void Geoave::push(double time,double e)
{
  int    n;
  double Q,R;

  if (time<t0) return;
  if (time==t0) n=0;
  else {
    errno=0;
    n= wfactor == 1 ? floor( (time-t0)/base) :
      floor( log(read_fb*(time-t0)+1)/logwf );
    if (errno) {
      std::cerr << "Time " << time << "(t0 = " <<t0<< ")\n";
      throw Runtime_error("Error in logarithm",HERE);
    }
    n+=1;
  }
  
  if (n>=count.size()) {
    count.resize(n+10,0);
    rave.resize(n+10,0);
    rvarn.resize(n+10,0);
  }

  count[n]++;
  Q=e-rave[n];
  R=Q/count[n];
  rave[n]+=R;
  rvarn[n]+=Q*R*(count[n]-1);
}

void Geoave::get_aves(std::vector<double>& time,std::vector<double>& ave,
		      std::vector<double>& var)
{
  time.clear();
  ave.clear();
  var.clear();

  double v,da;

  if (count[0]>0) {
    time.push_back(t0);
    ave.push_back(rave[0]);
    v=rvarn[0]/(count[0]-1);
    da=sqrt(v/count[0]);
    var.push_back(v);
  }
  
  double t=t0+0.5*base;
  double deltat=0.5*base*(1+wfactor);
  for (int i=1; i<rave.size(); i++, t+=deltat, deltat*=wfactor) {
    if (count[i]==0) continue;
    time.push_back(t);
    ave.push_back(rave[i]);
    v=rvarn[i]/(count[i]-1);
    da=sqrt(v/count[i]);
    var.push_back(v);
  }
}

std::ostream& operator<<(std::ostream& os,const Geoave& g)
{
  os << "# time   ave  deltaave(=s.d/sqrt(n))\n";
  double ave,var;
  if (g.count[0]>0) {
    var=g.rvarn[0]/(g.count[0]-1);
    os << g.t0 << "   " << g.rave[0] << "   " << sqrt(var/g.count[0]) << '\n';
  }

  double t=g.t0+0.5*g.base;
  double deltat=0.5*g.base*(1+g.wfactor);
  for (int i=1; i<g.rave.size(); i++, t+=deltat, deltat*=g.wfactor) {
    if (g.count[i]==0) continue;
    var=g.rvarn[i]/(g.count[i]-1);
    os << t << "  " << g.rave[i] << "  " << sqrt(var/g.count[i]) << '\n';
  }
  return os;
}


}
