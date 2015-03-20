/*
 * lattice1D.hh -- Lattices in 1 dimension
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

#ifndef _LATTICE1D_HH_
#define _LATTICE1D_HH_

#include <graph.hh>

/* @T

\section{Periodic 1-$d$ lattice}

A trivial 1-$d$ Bravais lattice with periodic boundary conditions,
mostly for pedagogical and testing purposes

*/

template <typename nodeT>
class Periodic1DLattice : public GraphBase<nodeT> {
public:
  Periodic1DLattice(int n);

  class node_iterator;
  class neighbour_iterator;

  node_iterator      begin() {return node_iterator(*this);}
  nodeT*             end() {return GraphBase<nodeT>::nodes+GraphBase<nodeT>::Nnodes;}
  neighbour_iterator neighbour_begin(ptrdiff_t id)
  {return neighbour_iterator(*this,id);}
  neighbour_iterator neighbour_begin(nodeT* n)
  {return neighbour_begin(this->id(n));}
  int neighbour_end(ptrdiff_t id) {return 0;}
  int neighbour_end(nodeT* n)  {return 0;}

} ;

template <typename nodeT> inline
Periodic1DLattice<nodeT>::Periodic1DLattice(int n) :
  GraphBase<nodeT>(n,2,false)
{}

// Node iterator
template <typename nodeT>
class Periodic1DLattice<nodeT>::node_iterator :
  public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:
  typedef Periodic1DLattice<nodeT> Graph_t;

  node_iterator(Periodic1DLattice<nodeT> &l) :
    lat(l), n(l.data()), LL(l.data()+l.size()-1), RR(l.data()+1) {}

  node_iterator(const node_iterator &i) :
  lat(i.lat), n(i.n), LL(i.LL), RR(i.RR)
  {}
  
  node_iterator& operator=(const node_iterator &it)
  {if (this==&it) return *this; lat=it.lat; n=it.n; RR=it.RR; LL=it.RR; return *this;}
  
  bool operator==(const node_iterator &i)
  {return n==i.n;}
  
  bool operator!=(const node_iterator &i)
  {return n!=i.n;}
  
  bool operator==(nodeT* p)
  {return n==p;}
  
  bool operator!=(nodeT* p)
  {return n!=p;}
  
  nodeT* L() const {return RR;}
  nodeT* R() const {return LL;}
  
  nodeT& operator*() const {return *n;}
  nodeT* operator->() const {return n;}
  
  node_iterator& operator++()
  {
    ++n; ++LL; ++RR;
    RR-=lat.size() * (lat.id(RR) / lat.size());
    return *this;
  }
  
  node_iterator operator++(int) {node_iterator c=*this; ++(*this); return c;}
  
  node_iterator& operator--()
  {
    --n;
    --LL;
    --RR;
    LL+=lat.size() * (lat.id(RR) / lat.size());
    return *this;
  }

  node_iterator operator--(int) {node_iterator c=*this; --(*this); return c;}

  // Above methods are required by a standard iterator; we add the following

  // node_iterator& to(int n) {nn=graph.nodes+nn; return *this;}
  // node_iterator& to_neighbour(int n) {nn=graph.node_neighbours(nn)+n; return *this;}
  // nodeT* neighbour(int i) const {return graph.node_neighbours(nn)+i;}

  // protected:
public:
  Graph_t&  lat;
  nodeT     *n,*LL,*RR;
} ;


// template <typename nodeT,typename Function>
// inline Function for_each_neighbour(PeriodicSQLattice<nodeT>& lat,
// 				   typename PeriodicSQLattice<nodeT>::node_iterator& n,Function f)
// {
//   f(n.N());
//   f(n.S());
//   f(n.E());
//   f(n.W());
// }



// Neighbour iterator

template <typename nodeT>
class Periodic1DLattice<nodeT>::neighbour_iterator :
  public GraphBase<nodeT>::neighbour_iterator
{
  using GraphBase<nodeT>::neighbour_iterator::nn;

public:
  neighbour_iterator(const Periodic1DLattice<nodeT>& graph,ptrdiff_t inode) :
  GraphBase<nodeT>::neighbour_iterator(nlist)
  {
    nlist[0]=inode==0 ? graph.nodes+graph.size()-1 : graph.nodes+inode-1;
    nlist[1]=inode==graph.size()-1 ? graph.nodes : graph.nodes+inode+1;
  }

  neighbour_iterator(const neighbour_iterator &ni) :
  GraphBase<nodeT>::neighbour_iterator(nlist)
  {nlist[0]=ni.nlist[0]; nlist[1]=ni.nlist[1];
    nn=nlist+(ni.nn-ni.nlist);}

  neighbour_iterator& operator=(const neighbour_iterator& ni)
  {nlist[0]=ni.nlist[0]; nlist[1]=ni.nlist[1];
    nn=nlist+(ni.nn-ni.nlist);}

  bool operator==(const typename GraphBase<nodeT>::neighbour_iterator &i)
  {return nn==i.nn;}

  bool operator!=(int p)
  {return nn!=nlist+2;}

  nodeT& operator*() const {return **nn;}
  nodeT* operator->() const {return *nn;}

  neighbour_iterator& operator++() {++nn; return *this;}
  neighbour_iterator& operator++(int) {neighbour_iterator c=*this; ++nn; return *this;}
  neighbour_iterator& operator--() {nn--; return *this;}
  neighbour_iterator& operator--(int) {neighbour_iterator c=*this; --nn; return *this;}

private:
  nodeT *nlist[2];
} ;

template <typename nodeT,typename Function>
inline Function for_each_neighbour(Periodic1DLattice<nodeT>::node_iterator& n,Function f)
{
  f(n.L());
  f(n.R());
}


#endif /* _LATTICE1D_HH_ */

