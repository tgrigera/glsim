/*
 * lattice2D.hh -- Lattices in two dimensions
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

#ifndef _LATTICE2D_HH_
#define _LATTICE2D_HH_

#include <math.h>

#include "graph.hh"

namespace glsim {

/** \ingroup lattice

 \brief The periodic square lattice (2 dimensions)

 This class implements a periodic SQ lattice.  The size \f$L_x \times
 L_y\f$ is specified on construction.  Sites can be accessed by a pair
 \f$(x,y)\f$ of coordinates, where \a x is the horizontal coordinate
 (if you prefer to see this as a matrix, then the first coordinate is
 the column number).

 Internally, data is stored by column (\a y is the fastest changing
 coordinate): this is the order you get from the `data()` method, or
 with the increment and decrement operations of the node_iterator.

 First neighbours are called north, south, east and west, and can be
 accessed by number or using these names with the corresponding
 node_iterator.

 The node_iterator includes the increment and decrement operators (as
 it should, since it is a bidirectional iterator.  However, these do
 not behave periodically but rather walk through all nodes and become
 equal to begin() or end().

*/
template <typename nodeT>
class PeriodicSQLattice : public GraphBase<nodeT> {
public:
  /// \name Construction, copy and destruction
  ///@{
  PeriodicSQLattice(int Lx,int Ly);
  PeriodicSQLattice(const PeriodicSQLattice &);
  ~PeriodicSQLattice();
  PeriodicSQLattice& operator=(const PeriodicSQLattice&);

  ///@}
  /// \name Lattice information
  ///@{

  int size_x() const {return Lx;}
  int size_y() const {return Ly;}
  id_t         id(int x,int y) const {return lpos(x,y);}
  id_t         id_protected(int x,int y) const {
                                         return lpos(mymod(x,Lx),mymod(y,Ly));}
  ///@}
  /// \name Node access
  ///@{
  using GraphBase<nodeT>::id;
  using GraphBase<nodeT>::data;
  nodeT&       operator()(id_t x,id_t y) {return (data())[lpos(x,y)];}
  const nodeT& operator()(id_t x,id_t y) const {return (data())[lpos(x,y)];}

  ///@}
  ///\name Disk i/o
  ///@{
  void read(std::istream&);  ///< Read graph from stream
  void write(std::ostream&); ///< Write graph to stream

  ///@}

  ///@}
  ///\name Utility
  ///@{
  int pdiff(int i,int j,int L) const;  ///<Periodic difference
  int pdiffx(int i, int j) const {return pdiff(i,j,Lx);}
  int pdiffy(int i, int j) const {return pdiff(i,j,Ly);}

  ///@}
  /// \name Iterators
  ///@{
  class   node_iterator;

  node_iterator      begin() {return node_iterator(*this);}
  nodeT*             end()   {return GraphBase<nodeT>::end();}
  ///@}

private:
  int  Lx,Ly;
  int  *Ndis,*Sdis,*Edis,*Wdis;   // To implement periodic conditions

  // This transforms coordinates to a valid node id
  id_t lpos(id_t i,id_t j) const {return Ly*i + j;}
  void posl(id_t l,id_t &i,id_t &j) const;

  int mymod(int a,int b) {int r=a%b; return (r>=0) ? r : r+b;}
  void load_neighbour_distance_vectors();
} ;

/*******************************************************************************
 * Construction, copy and destruction
 *
 */

/** \brief Create a square lattice of size \f$L_x \times L_y\f$
 */
template <typename nodeT>
inline PeriodicSQLattice<nodeT>::PeriodicSQLattice(int Lx_,int Ly_) :
  GraphBase<nodeT>(Lx_*Ly_,4,false),
  Lx(Lx_), Ly(Ly_),
  Ndis(0), Sdis(0), Edis(0), Wdis(0)
{
  load_neighbour_distance_vectors();
}

template <typename nodeT>
inline PeriodicSQLattice<nodeT>::PeriodicSQLattice(const PeriodicSQLattice &l) :
  GraphBase<nodeT>(l),
  Lx(l.Lx), Ly(l.Ly),
  Ndis(0), Sdis(0), Edis(0), Wdis(0)
{
  load_neighbour_distance_vectors();
}

template<typename nodeT>
PeriodicSQLattice<nodeT>& PeriodicSQLattice<nodeT>::
operator=(const PeriodicSQLattice<nodeT>& l)
{
  if (this==&l) return *this;
  GraphBase<nodeT>::operator=(l);
  Lx=l.Lx;
  Ly=l.Ly;
  load_neighbour_distance_vectors();
  return *this;
}

template <typename nodeT>
void PeriodicSQLattice<nodeT>::load_neighbour_distance_vectors()
{
  int i,j;

  if (Ndis) delete[] Ndis;
  if (Sdis) delete[] Sdis;
  if (Edis) delete[] Edis;
  if (Wdis) delete[] Wdis;
  Ndis=new int[Ly];
  Sdis=new int[Ly];
  Wdis=new int[Lx];
  Edis=new int[Lx];

  for (i=0; i<Ly; i++) {
    Ndis[i]=1;
    Sdis[i]=-1;
  }
  Sdis[0]=Ly-1;
  Ndis[Ly-1]=-Sdis[0];
  for (i=0; i<Lx; i++) {
    Edis[i]=Ly;
    Wdis[i]=-Ly;
  }
  Wdis[0]=Ly*(Lx-1);
  Edis[Lx-1]=-Wdis[0];
}

template<typename nodeT>
PeriodicSQLattice<nodeT>::~PeriodicSQLattice()
{
  delete[] Ndis;
  delete[] Sdis;
  delete[] Wdis;
  delete[] Edis;
}

/******************************************************************************
 *
 * Utilitty and i/o
 *
 */

template<typename nodeT> inline
void PeriodicSQLattice<nodeT>::posl(id_t l,id_t &i,id_t &j) const
{
  i=l/Ly;
  j=l-i*Ly;
}

template<typename nodeT>
inline int PeriodicSQLattice<nodeT>::pdiff(int i,int j,int L) const
{
  return i-j-L*rint((double)(i-j)/L);
}

template<typename nodeT>
void PeriodicSQLattice<nodeT>::write(std::ostream& f)
{
  GraphBase<nodeT>::write(f);
  f.write((char*) &Lx,sizeof(Lx));
  f.write((char*) &Ly,sizeof(Ly));
}

template<typename nodeT>
void PeriodicSQLattice<nodeT>::read(std::istream& f)
{
  GraphBase<nodeT>::read(f);
  f.read((char*) &Lx,sizeof(Lx));
  f.read((char*) &Ly,sizeof(Ly));
  load_neighbour_distance_vectors();
}

/******************************************************************************
 *
 * Iterator
 *
 */
/** \brief A bidirectional iterator for the periodic SQ lattice

 This is a neighbour-aware iterator that works like
 GraphBase<nodeT>::node_iterator, except that it adds the methods N(),
 S(), E() and W() to access the first neighbours of a node (through a
 reference).  For site \f$(i,j)\f$, the North neighbour is
 \f$(i,j+1)\f$, East is \f$(i+1,j)\f$, etc.
 
                                N
                           W  (i,j)  E
                                S

 Neighbours can also be accessed by number (from 0 to 3), see neighbour().

 To move the iterator (make it point to another node), appart from the
 standard bidirectional iterator methods you can use to() (overloaded
 to accept coordinates or a pointer) and to_neighbour().

 */
template <typename nodeT>
class PeriodicSQLattice<nodeT>::node_iterator :
  public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:
  typedef PeriodicSQLattice<nodeT> Graph_t;  ///< The type of the graph we belong to

  ///\name Construction, copy and comparison
  //@{
  node_iterator(PeriodicSQLattice<nodeT> &l,id_t i_=0,id_t j_=0);
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


/*******************************************************************************
 *
 * Construction, destruction and comparison
 *
 */

template <typename nodeT> inline PeriodicSQLattice<nodeT>::node_iterator::
node_iterator(PeriodicSQLattice<nodeT> &lattice ,id_t i_,id_t j_) :
 lat(lattice),
 i(i_),
 j(j_),
 n(lattice.data()+lattice.id(i_,j_))
{}

template <typename nodeT> inline PeriodicSQLattice<nodeT>::node_iterator::
node_iterator(PeriodicSQLattice<nodeT> &lattice,nodeT* np) :
  lat(lattice),
  n(np)
{
  lat.posl(lat.id(n),i,j);
}

template <typename nodeT> inline PeriodicSQLattice<nodeT>::node_iterator::
node_iterator(const node_iterator &i) :
 lat(i.lat),
 i(i.i),
 j(i.j),
 n(i.n)
{}

template <typename nodeT> inline
typename PeriodicSQLattice<nodeT>::node_iterator&
PeriodicSQLattice<nodeT>::node_iterator::operator=(const node_iterator &it)
{
  if (this==&it) return *this;
  lat=it.lat;
  i=it.i;
  j=it.j; 
  n=it.n;
  return *this;
}

/*
 * Movement and neighbours
 *
 */

template <typename nodeT> inline
typename PeriodicSQLattice<nodeT>::node_iterator&
PeriodicSQLattice<nodeT>::node_iterator::to(id_t i_,id_t j_)
{
  i=i_;
  j=j_;
  n=lat.data()+lat.id(i,j);
  return *this;
}

template <typename nodeT>
typename PeriodicSQLattice<nodeT>::node_iterator&
PeriodicSQLattice<nodeT>::node_iterator::to(nodeT* node)
{
  id_t i,j;
  lat.posl(lat.id(node),i,j);
  return to(i,j);
}
 
template <typename nodeT>
nodeT& PeriodicSQLattice<nodeT>::node_iterator::neighbour(int i) const
{
  switch(i) {
  case 0: return N();
  case 1: return E();
  case 2: return S();
  case 3: return W();
  }
}

template <typename nodeT> inline
typename PeriodicSQLattice<nodeT>::node_iterator&
PeriodicSQLattice<nodeT>::node_iterator::to_neighbour(int i)
{
  switch(i) {
  case 0: return to(&N());
  case 1: return to(&E());
  case 2: return to(&S());
  case 3: return to(&W());
  }
  throw glsim::Internal_error(HERE);
}


template <typename nodeT,typename Function>
struct implement_for_each_neighbour<nodeT,Function,typename PeriodicSQLattice<nodeT>::node_iterator> {
  inline static Function fen(typename PeriodicSQLattice<nodeT>::node_iterator n,Function f)
  {
    f(n.N());
    f(n.E());
    f(n.S());
    f(n.W());
    return f;
  }
} ;


template <typename nodeT,typename Function>
struct implement_for_each_bond<nodeT,Function,PeriodicSQLattice<nodeT>> {
  inline static Function fen(PeriodicSQLattice<nodeT> &lat,Function f)
  {
    for (auto iter=lat.begin(); iter!=lat.end(); ++iter) {
      f(*iter,iter.W());
      f(*iter,iter.S());
    }
    return f;
  }
} ;

} /* namespace */

#endif /* _LATTICE2D_HH_ */
