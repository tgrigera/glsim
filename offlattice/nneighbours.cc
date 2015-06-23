/*
 * nnneighbours.cc -- structures to find nearest neighbours
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

#include "nneighbours.hh"

namespace glsim {

/*****************************************************************************
 *
 * NeighbourList_naive
 *
 */


NeighbourList_naive::NeighbourList_naive(double rc_,double delta_r_) :
  NearestNeighbours(rc_),
  topological(false), Nnearest(0),
  delta_r(delta_r_),
  accum_maxdisp(0.)
{
  if (delta_r<0)
    delta_r=rc*0.3;
  rdsq=rc+delta_r;
  rdsq*=rdsq;
}

void NeighbourList_naive::rebuild(OLconfiguration& conf,double rc_)
{
  if (rc_>0) {
    rc=rc_;
    rdsq=rc+delta_r;
    rdsq*=rdsq;
  }
  topological=false;
  
  neighbours.clear();
  neighbours.resize(conf.N);
  pairs.clear();
  for (size_t i=0; i<conf.N-1; ++i)
    for (size_t j=i+1; j<conf.N; ++j)
      if (conf.distancesq(i,j) <= rdsq) {
	pairs.push_back(Pair(i,j));
	neighbours[i].push_back(j);
	neighbours[j].push_back(i);
      }
  accum_maxdisp=0.;
}

void NeighbourList_naive::rebuild(OLconfiguration& conf,int n)
{
  Nnearest=n;
  topological=true;
  
  throw glsim::Unimplemented("Topological neighbours");
}

void NeighbourList_naive::update(OLconfiguration& conf,double maxdisp)
{
  if (topological)
    rebuild(conf,Nnearest);
  else {
    accum_maxdisp+=maxdisp;
    if (accum_maxdisp>=delta_r/2.)
      rebuild(conf);
  }
}

} /* namespace */
