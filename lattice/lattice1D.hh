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

#include "graph.hh"

namespace glsim {

/** \ingroup lattice

 \brief  A periodic lattice in one dimension.

 A trivial 1-d Bravais lattice with periodic boundary conditions,
 mostly for pedagogical and testing purposes.  Provides a
 node_iterator that knows about left and right neighbours.

*/
template <typename nodeT>
class Periodic1DLattice : public GraphBase<nodeT> {
public:
  Periodic1DLattice(int N);

  class node_iterator;

  node_iterator begin() {return node_iterator(*this);}
  nodeT*        end() {return GraphBase<nodeT>::data()+GraphBase<nodeT>::size();}
} ;

/** \brief Create a periodic lattice with the specified number of nodes
 */
template <typename nodeT> inline
Periodic1DLattice<nodeT>::Periodic1DLattice(int N) :
  GraphBase<nodeT>(N,2,false)
{}

/** \brief Bidirectional iterator for Periodic1DLattice
 */
template <typename nodeT>
class Periodic1DLattice<nodeT>::node_iterator :
  public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:
  typedef Periodic1DLattice<nodeT> Graph_t; ///< The type of the graph we belong to

  ///\name Construction, copy and comparison
  //@{
  node_iterator(Periodic1DLattice<nodeT> &l,id_t i=0);
  node_iterator(Periodic1DLattice<nodeT> &l,nodeT*);

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
  {
    ++n;
    ++(ng[0]);
    ++(ng[1]);
    ng[1]-=lat.size() * (lat.id(ng[1]) / lat.size());
    return *this;
  }
  
  node_iterator operator++(int) {node_iterator c=*this; ++(*this); return c;}
  
  node_iterator& operator--()
  {
    --n;
    --(ng[0]);
    --(ng[1]);
    ng[0]+=lat.size() * (lat.id(ng[0]) / lat.size());
    return *this;
  }

  node_iterator operator--(int) {node_iterator c=*this; --(*this); return c;}

  ///@}

  ///\name Extra methods (including neighbour access)
  ///@{
  node_iterator& to(id_t i);
  node_iterator& to(nodeT *n_) {to(lat.id(n_));}

  int            neighbour_size() const {return 2;}
  node_iterator& to_neighbour(int n);

  nodeT&         neighbour(int i) const {return *(ng[i]);}
  nodeT&         L() const {return *(ng[0]);}
  nodeT&         R() const {return *(ng[1]);}

  operator nodeT*() const {return n;}

  ///@}

private:
  Graph_t &lat;
  nodeT   *n,*ng[2];
} ;


template <typename nodeT> inline
Periodic1DLattice<nodeT>::node_iterator::
node_iterator(Periodic1DLattice<nodeT> &l,id_t i) :
  lat(l)
{
  to(i);
}

template <typename nodeT> inline
Periodic1DLattice<nodeT>::node_iterator::
node_iterator(Periodic1DLattice<nodeT> &l,nodeT *n_) :
  lat(l)
{
  to(n_);
}

template <typename nodeT> inline
Periodic1DLattice<nodeT>::node_iterator::
node_iterator(const node_iterator &i) :
  lat(i.lat),
  n(i.n)
{
  ng[0]=i.ng[0];
  ng[1]=i.ng[1];
}

template <typename nodeT> inline
typename Periodic1DLattice<nodeT>::node_iterator&
Periodic1DLattice<nodeT>::node_iterator::
operator=(const node_iterator &it)
{
  if (this==&it) return *this;
  lat=it.lat;
  n=it.n;
  ng[0]=it.ng[0];
  ng[1]=it.ng[1];
  return *this;
}

/*
 * movement
 *
 */

template <typename nodeT>
typename Periodic1DLattice<nodeT>::node_iterator&
Periodic1DLattice<nodeT>::node_iterator::to(id_t i)
{
  n=lat.data()+i;
  ng[0]=n-1;
  ng[0]+=lat.size() * (lat.id(ng[0]) / lat.size());
  ng[1]=n+1;
  ng[1]-=lat.size() * (lat.id(ng[1]) / lat.size());
  return *this;
}

template <typename nodeT> inline
typename Periodic1DLattice<nodeT>::node_iterator&
Periodic1DLattice<nodeT>::node_iterator::to_neighbour(int i)
{
  id_t nid=GraphBase<nodeT>::id(ng[i]);
  to(nid);
}

// template <typename nodeT,typename Function>
// inline Function for_each_neighbour(typename Periodic1DLattice<nodeT>::node_iterator& n,Function f)
// {
// std::cout << "This is meee\n"; exit(0);
//   f(n.L());
//   f(n.R());
// }


} /* namespace */

#endif /* _LATTICE1D_HH_ */

