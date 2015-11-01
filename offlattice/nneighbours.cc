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

void NeighbourList_naive::update(OLconfiguration& conf,double maxdisp)
{
  accum_maxdisp+=maxdisp;
  if (accum_maxdisp>=delta_r/2.)
    rebuild(conf);
}

/*****************************************************************************
 *
 * TopologialNeighbours_naive
 *
 */

void TopologicalNeighbours_naive::rebuild(OLconfiguration& conf,int Nnearest_)
{
  if (Nnearest_>0) Nnearest=Nnearest_;
}

vu()
{
  neighbours.clear();
  neighbours.resize(conf.N);
  pairs.clear();


  std::priority_queue<ndist,std::vector<ndist>,comp> candidates;

  for (size_t i=0; i<conf.N; ++i) {
    candidates.clear();
    for (size_t j=0; j<conf.N; ++j) {
      if (i==j) continue;
      double dsq=conf.distancesq(i,j);
      if (candidates.size()<Nnearest)
	candidates.push(ndist(j,dsq));
      else if (dsq<candidates.top().dsq) {
	candidates.pop();
	candidates.push(ndist(j,dsq));
      }
    }
  }
	

// 	pairs.push_back(Pair(i,j));
// 	neighbours[i].push_back(j);
// 	neighbours[j].push_back(i);
//       }
// }

}

#ifdef UFA
class TopologicalNeighbours_naive::neighbour_iterator :
  public std::iterator<std::forward_iterator_tag,int>
{
public:
  ///\name Construction, copy and comparison
  //@{
  nieghbour_iterator(PeriodicSQLattice<nodeT> &l,id_t i_=0,id_t j_=0);
  ///< Create and move to node (i,j)
  node_iterator(PeriodicSQLattice<nodeT> &l,nodeT* n);
  ///< Create and move to node n
  node_iterator(const node_iterator &i);

  node_iterator& operator=(const node_iterator &it);

  bool operator==(const node_iterator &i)
  {return n==i.n;}

  bool operator!=(const node_iterator &i)
  {return n!=i.n;}

  bool operator==(nodeT* p)
  {return n==p;}

  bool operator!=(nodeT* p)
  {return n!=p;}

  //@}
  ///\name Operators required by standard bidirectional iterators
  ///@{
  nodeT& operator*() const {return *n;}

  nodeT* operator->() const {return n;}

  node_iterator& operator++()
  {++n; i+=j/(lat.Ly-1); j+=lat.Ndis[j]; return *this;}

  node_iterator operator++(int) {node_iterator c=*this; ++(*this); return c;}

  node_iterator& operator--()
  {--n; j+=lat.Sdis[j]; i-=j/(lat.Lx-1);   return *this;}

  node_iterator operator--(int) {node_iterator c=*this; --(*this); return c;}
  //@}

  ///\name Extra methods (including neighbour access)
  ///@{
  node_iterator& to(id_t i,id_t j);  ///< Move to node (i,j)
  node_iterator& to(nodeT* n);       ///< Move to node specified by pointer

  int            neighbour_size() const {return 4;}  ///<Coordination number of node (always 4)
  node_iterator& to_neighbour(int n);  ///< Move to nth neighbour (numbered as in neighbour())

  nodeT& neighbour(int i) const; ///< Access neighbours by number (0=north, then clockwise up to 3=west)
  nodeT& N() const {return *(n+lat.Ndis[j]);}
  nodeT& S() const {return *(n+lat.Sdis[j]);}
  nodeT& E() const {return *(n+lat.Edis[i]);}
  nodeT& W() const {return *(n+lat.Wdis[i]);}

  id_t  x() {return i;}
  id_t  y() {return j;}
  operator nodeT*() const {return n;}

  ///@}

private:
  Graph_t&  lat;
  id_t      i,j;
  nodeT*    n;
} ;

#endif













} /* namespace */
