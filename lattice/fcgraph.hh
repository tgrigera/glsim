/*
 * fcgraph.hh -- classes for the Bethe lattice
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019, 2020, 2021
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

#ifndef FCGRAPH_HH
#define FCGRAPH_HH

#include "graph.hh"

namespace glsim {

/** \ingroup lattice
 * \brief class for the fully-connected graph
 *
 * The fully-connected graph has edges connecting all pairs of nodes
 *
 */

template <typename nodeT>
class FullyConnectedGraph : public GraphBase<nodeT> {
public:
  FullyConnectedGraph(int N); 

  ///@}
  /// \name Iterators
  ///@{
  class   node_iterator;

  node_iterator      begin() {return node_iterator(*this);}
  nodeT*             end()   {return GraphBase<nodeT>::end();}
  ///@}

private:

} ;
  
template <typename nodeT>
FullyConnectedGraph<nodeT>::FullyConnectedGraph(int N) :
  GraphBase<nodeT>(N,N-1,false)
{
}


/******************************************************************************
 *
 * Iterator
 *
 */
/** \brief A bidirectional iterator for the fully-connected graph

 This is a neighbour-aware iterator that works like
 GraphBase<nodeT>::node_iterator.

 Neighbours can also be accessed by number (from 0 to N-1), see neighbour().

 To move the iterator (make it point to another node), appart from the
 standard bidirectional iterator methods you can use to() (overloaded
 to accept coordinates or a pointer) and to_neighbour().

 */
template <typename nodeT>
class FullyConnectedGraph<nodeT>::node_iterator :
  public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:
  typedef FullyConnectedGraph<nodeT> Graph_t;  ///< The type of the graph we belong to

  ///\name Construction, copy and comparison
  //@{
  node_iterator(FullyConnectedGraph<nodeT> &l,id_t n_=0);
  ///< Create and move to node (i,j)
  node_iterator(FullyConnectedGraph<nodeT> &l,nodeT* n);
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
  {++n; return *this;}

  node_iterator operator++(int) {node_iterator c=*this; ++(*this); return c;}

  node_iterator& operator--()
  {--n; return *this;}

  node_iterator operator--(int) {node_iterator c=*this; --(*this); return c;}
  //@}

  ///\name Extra methods (including neighbour access)
  ///@{
  node_iterator& to(id_t n);  ///< Move to node n
  node_iterator& to(nodeT* n);       ///< Move to node specified by pointer

  int            neighbour_size() const {return lat.coordination_number();}  ///<Coordination number of node
  node_iterator& to_neighbour(int n);  ///< Move to nth neighbour (numbered as in neighbour())

  nodeT& neighbour(int i) const; ///< Access neighbours by number

  operator nodeT*() const {return n;}

  ///@}

private:
  Graph_t&  lat;
  nodeT*    n;
} ;

/*******************************************************************************
 *
 * Construction, destruction and comparison
 *
 */

template <typename nodeT> inline FullyConnectedGraph<nodeT>::node_iterator::
node_iterator(FullyConnectedGraph<nodeT> &lattice ,id_t n_) :
 lat(lattice),
 n(lattice.data()+n_)
{}

template <typename nodeT> inline FullyConnectedGraph<nodeT>::node_iterator::
node_iterator(FullyConnectedGraph<nodeT> &lattice,nodeT* np) :
  lat(lattice),
  n(np)
{}

template <typename nodeT> inline FullyConnectedGraph<nodeT>::node_iterator::
node_iterator(const node_iterator &i) :
 lat(i.lat),
 n(i.n)
{}

template <typename nodeT> inline
typename FullyConnectedGraph<nodeT>::node_iterator&
FullyConnectedGraph<nodeT>::node_iterator::operator=(const node_iterator &it)
{
  if (this==&it) return *this;
  lat=it.lat;
  n=it.n;
  return *this;
}

/*
 * Movement and neighbours
 *
 */

template <typename nodeT> inline
typename FullyConnectedGraph<nodeT>::node_iterator&
FullyConnectedGraph<nodeT>::node_iterator::to(id_t n_)
{
  n=lat.data()+n_;
  return *this;
}

template <typename nodeT>
typename FullyConnectedGraph<nodeT>::node_iterator&
FullyConnectedGraph<nodeT>::node_iterator::to(nodeT* node)
{
  n=node;
  return *this;
}
 
template <typename nodeT>
nodeT& FullyConnectedGraph<nodeT>::node_iterator::neighbour(int i) const
{
  nodeT *neigh=lat.data()+i;
  if (neigh>=n) ++neigh;
  return *neigh;
}

template <typename nodeT> inline
typename FullyConnectedGraph<nodeT>::node_iterator&
FullyConnectedGraph<nodeT>::node_iterator::to_neighbour(int i)
{
  return to(&neighbour(i));
}

template <typename nodeT,typename Function>
struct implement_for_each_neighbour<nodeT,Function,typename FullyConnectedGraph<nodeT>::node_iterator> {
  inline static Function fen(typename FullyConnectedGraph<nodeT>::node_iterator n,Function f)
  {
    for (int i=0; i<n.neighbour_size(); ++i)
      f(n.neighbour(i));
    return f;
  }
} ;


template <typename nodeT,typename Function>
struct implement_for_each_bond<nodeT,Function,FullyConnectedGraph<nodeT>> {
  inline static Function fen(FullyConnectedGraph<nodeT> &lat,Function f)
  {
    nodeT *end=lat.data()+lat.size();
    for (nodeT *n1=lat.data(); n1!=end-1; ++n1)
      for (nodeT *n2=n1+1; n2!=end; ++n2)
	f(*n1,*n2);
    return f;
  }
} ;

} /* namespace */

#endif /* FCGRAPH_HH */
