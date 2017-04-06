/*
 * lattice3D.hh -- Lattices in three dimensions
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

#ifndef LATTICE3D_HH
#define LATTICE3D_HH

#include <math.h>

#include "graph.hh"

namespace glsim {

/** \ingroup lattice

 \brief The periodic simple cubic lattice (3 dimensions).

 This class implements a simple cubic (SC) lattice with periodic
 boundary conditions.  The size \f$L_x \times L_y \times L_z\f$ is
 specified on construction.  Sites can be accessed by specifiying
 three coordinates \f$(x,y,z)\f$.

 Internally, data is stored so that z is the fastest changing
 coordinate: this is the order you get from the `data()` method, or
 with the increment and decrement operations of the node_iterator.

 First neighbours are called north, south, east, west, up, and down
 and can be accessed by number or using these names with the
 corresponding node_iterator.

 The node_iterator includes the increment and decrement operators (as
 it should, since it is a bidirectional iterator.  However, these do
 not behave periodically but rather walk through all nodes and become
 equal to begin() or end().

*/
template <typename nodeT>
class PeriodicSCLattice : public GraphBase<nodeT> {
public:
  /// \name Construction, copy and destruction
  ///@{
  PeriodicSCLattice(int Lx,int Ly,int Lz);
  PeriodicSCLattice(const PeriodicSCLattice &);
  ~PeriodicSCLattice();
  PeriodicSCLattice& operator=(const PeriodicSCLattice&);

  ///@}
  /// \name Lattice information
  ///@{
  int size_x() const {return Lx;}
  int size_y() const {return Ly;}
  int size_z() const {return Lz;}
  id_t         id(int x,int y,int z) const {return lpos(x,y,z);}
  id_t         id_protected(int x,int y,int z) const
  {return lpos(mymod(x,Lx),mymod(y,Ly));}

  ///@}
  /// \name Node access
  ///@{
  using GraphBase<nodeT>::id;
  using GraphBase<nodeT>::data;
  nodeT&       operator()(int x,int y,int z) {return (data())[lpos(x,y,z)];}
  const nodeT& operator()(int x,int y,int z) const {return (data())[lpos(x,y,z)];}

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
  int pdiffz(int i, int j) const {return pdiff(i,j,Lz);}

  ///@}
  /// \name Iterators
  ///@{
  class   node_iterator;

  node_iterator      begin() {return node_iterator(*this);}
  nodeT*             end()   {return GraphBase<nodeT>::end();}
  ///@}

private:
  int  Lx,Ly,Lz;
  int  *Ndis,*Sdis,*Edis,*Wdis,*Udis,*Ddis;   // To implement periodic conditions

  // This transforms coordinates to a valid node id
  id_t lpos(int i,int j,int k) const {return Lz*Ly*i + Lz*j + k;}
  void posl(id_t l,int &i,int &j,int &k) const;

  int mymod(int a,int b) {int r=a%b; return (r>=0) ? r : r+b;}
  void load_neighbour_distance_vectors();
} ;

/*******************************************************************************
 * Construction, copy and destruction
 *
 */

/** \brief Create a simple cubic lattice of size \f$L_x \times L_y\times L_z\f$
 */
template <typename nodeT>
inline PeriodicSCLattice<nodeT>::PeriodicSCLattice(int Lx_,int Ly_,int Lz_) :
  GraphBase<nodeT>(Lx_*Ly_*Lz_,6,false),
  Lx(Lx_), Ly(Ly_), Lz(Lz_),
  Ndis(0), Sdis(0), Edis(0), Wdis(0), Udis(0), Ddis(0)
{
  load_neighbour_distance_vectors();
}

template <typename nodeT>
inline PeriodicSCLattice<nodeT>::PeriodicSCLattice(const PeriodicSCLattice &l) :
  GraphBase<nodeT>(l),
  Lx(l.Lx), Ly(l.Ly), Lz(l.Lz),
  Ndis(0), Sdis(0), Edis(0), Wdis(0), Udis(0), Ddis(0)
{
  load_neighbour_distance_vectors();
}

template<typename nodeT>
PeriodicSCLattice<nodeT>& PeriodicSCLattice<nodeT>::
operator=(const PeriodicSCLattice<nodeT>& l)
{
  if (this==&l) return *this;
  GraphBase<nodeT>::operator=(l);
  Lx=l.Lx;
  Ly=l.Ly;
  Lz=l.Lz;
  load_neighbour_distance_vectors();
  return *this;
}

template <typename nodeT>
void PeriodicSCLattice<nodeT>::load_neighbour_distance_vectors()
{
  int i,j;

  if (Ndis) delete[] Ndis;
  if (Sdis) delete[] Sdis;
  if (Edis) delete[] Edis;
  if (Wdis) delete[] Wdis;
  if (Ddis) delete[] Ddis;
  if (Udis) delete[] Udis;
  Ndis=new int[Ly];
  Sdis=new int[Ly];
  Wdis=new int[Lx];
  Edis=new int[Lx];
  Udis=new int[Lz];
  Ddis=new int[Lz];

  for (i=0; i<Lz; i++) {
    Udis[i]=1;
    Ddis[i]=-1;
  }
  Ddis[0]=Lz-1;
  Udis[Lz-1]=-Ddis[0];
  for (i=0; i<Ly; i++) {
    Ndis[i]=Lz;
    Sdis[i]=-Lz;
  }
  Sdis[0]=Lz*(Ly-1);
  Ndis[Ly-1]=-Sdis[0];
  for (i=0; i<Lx; i++) {
    Edis[i]=Ly*Lz;
    Wdis[i]=-Ly*Lz;
  }
  Wdis[0]=Lz*Ly*(Lx-1);
  Edis[Lx-1]=-Wdis[0];
}

template<typename nodeT>
PeriodicSCLattice<nodeT>::~PeriodicSCLattice()
{
  delete[] Ndis;
  delete[] Sdis;
  delete[] Wdis;
  delete[] Edis;
  delete[] Udis;
  delete[] Ddis;
}

/******************************************************************************
 *
 * Utilitty and i/o
 *
 */

template<typename nodeT> inline
void PeriodicSCLattice<nodeT>::posl(id_t l,int &i,int &j,int &k) const
{
  i=l/(Ly*Lz);
  j=(l-i*Ly*Lz)/Lz;
  k=l-i*Ly*Lz-j*Lz;
}

template<typename nodeT>
inline int PeriodicSCLattice<nodeT>::pdiff(int i,int j,int L) const
{
  return i-j-L*rint((double)(i-j)/L);
}

template<typename nodeT>
void PeriodicSCLattice<nodeT>::write(std::ostream& f)
{
  GraphBase<nodeT>::write(f);
  f.write((char*) &Lx,sizeof(Lx));
  f.write((char*) &Ly,sizeof(Ly));
  f.write((char*) &Lz,sizeof(Lz));
}

template<typename nodeT>
void PeriodicSCLattice<nodeT>::read(std::istream& f)
{
  GraphBase<nodeT>::read(f);
  f.read((char*) &Lx,sizeof(Lx));
  f.read((char*) &Ly,sizeof(Ly));
  f.read((char*) &Lz,sizeof(Lz));
  load_neighbour_distance_vectors();
}

/******************************************************************************
 *
 * Iterator
 *
 */
/** \brief A bidirectional iterator for the periodic SC lattice

 This is a neighbour-aware iterator that works like
 GraphBase<nodeT>::node_iterator, except that it adds the methods N(),
 S(), E(), W(), U(), and D() to access the first neighbours of a node
 (through a reference).  For site \f$(i,j,k)\f$, the North neighbour is
 \f$(i,j+1,k)\f$, East is \f$(i+1,j,k)\f$,  Up is \f$(i,j,k+1)\f$, etc.
 
                                 N
                           W  (i,j,k)  E
                                 S

 Neighbours can also be accessed by number (from 0 to 5), see neighbour().

 To move the iterator (make it point to another node), appart from the
 standard bidirectional iterator methods you can use to() (overloaded
 to accept coordinates or a pointer) and to_neighbour().

 */
template <typename nodeT>
class PeriodicSCLattice<nodeT>::node_iterator :
  public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:
  typedef PeriodicSCLattice<nodeT> Graph_t;  ///< The type of the graph we belong to

  ///\name Construction, copy and comparison
  //@{
  node_iterator(PeriodicSCLattice<nodeT> &l,int i_=0,int j_=0,int k_=0);
  ///< Create and move to (i,j,k)
  node_iterator(PeriodicSCLattice<nodeT> &l,nodeT* n);
  ///< Create and move to n

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
  {++n; i+=(j/(lat.Ly-1))*(k/(lat.Lz-1));
    j+=(k/(lat.Lz-1))*lat.Ndis[j]/lat.Lz; k+=lat.Udis[k]; return *this;}

  node_iterator operator++(int) {node_iterator c=*this; ++(*this); return c;}

  node_iterator& operator--()
  {--n; k+=lat.Ddis[k]; j+=(k/lat.Lz-1)*lat.Sdis[j]/lat.Lz;
    i-=(j/(lat.Ly-1))*(k/(lat.Lz-1)); return *this;}

  node_iterator operator--(int) {node_iterator c=*this; --(*this); return c;}
  //@}

  ///\name Extra methods (including neighbour access)
  ///@{
  node_iterator& to(int i,int j,int k); ///< Move to node given by coordinates.
  node_iterator& to(nodeT* n);  ///< Move to node given by pointer.

  int             neighbour_size() const {return 6;}
  ///< Coordination number of node (always 6).
  node_iterator& to_neighbour(int i); ///<Move to `i`th neighbour (numbers as in neighbour()).

  nodeT& neighbour(int i) const;
  ///< Access neighbour by number (0=north, clockwise to 3=west, 4=up, 5=down).
  nodeT& N() const {return *(n+lat.Ndis[j]);}
  nodeT& S() const {return *(n+lat.Sdis[j]);}
  nodeT& E() const {return *(n+lat.Edis[i]);}
  nodeT& W() const {return *(n+lat.Wdis[i]);}
  nodeT& U() const {return *(n+lat.Udis[k]);}
  nodeT& D() const {return *(n+lat.Ddis[k]);}

  int  x() {return i;}
  int  y() {return j;}
  int  z() {return k;}
  operator nodeT*() const {return n;}

  ///@}

private:
  Graph_t&  lat;
  int       i,j,k;
  nodeT*    n;
} ;


/*******************************************************************************
 *
 * Construction, destruction and comparison
 *
 */

template <typename nodeT> inline PeriodicSCLattice<nodeT>::node_iterator::
node_iterator(PeriodicSCLattice<nodeT> &lattice ,int i_,int j_,int k_) :
 lat(lattice),
 i(i_),
 j(j_),
 k(k_),
 n(lattice.data()+lattice.id(i_,j_,k_))
{}

template <typename nodeT> inline PeriodicSCLattice<nodeT>::node_iterator::
node_iterator(PeriodicSCLattice<nodeT> &lattice ,nodeT* n_) :
 lat(lattice),
 n(n_)
{
  lat.posl(lat.id(n),i,j,k);
}
  
template <typename nodeT> inline PeriodicSCLattice<nodeT>::node_iterator::
node_iterator(const node_iterator &i) :
 lat(i.lat),
 i(i.i),
 j(i.j),
 k(i.k),
 n(i.n)
{}

template <typename nodeT> inline
typename PeriodicSCLattice<nodeT>::node_iterator&
PeriodicSCLattice<nodeT>::node_iterator::operator=(const node_iterator &it)
{
  if (this==&it) return *this;
  lat=it.lat;
  i=it.i;
  j=it.j; 
  k=it.k;
  n=it.n;
  return *this;
}

/*
 * Movement and neighbours
 *
 */

template <typename nodeT> inline
typename PeriodicSCLattice<nodeT>::node_iterator&
PeriodicSCLattice<nodeT>::node_iterator::to(int i_,int j_,int k_)
{
  i=i_;
  j=j_;
  k=k_;
  n=lat.data()+lat.id(i,j,k);
  return *this;
}

template <typename nodeT>
typename PeriodicSCLattice<nodeT>::node_iterator&
PeriodicSCLattice<nodeT>::node_iterator::to(nodeT* node)
{
  int i,j,k;
  lat.posl(lat.id(node),i,j,k);
  return to(i,j,k);
}
 
template <typename nodeT>
nodeT& PeriodicSCLattice<nodeT>::node_iterator::neighbour(int i) const
{
  switch(i) {
  case 0: return N();
  case 1: return E();
  case 2: return S();
  case 3: return W();
  case 4: return U();
  case 5: return D();
  }
}

template <typename nodeT> inline
typename PeriodicSCLattice<nodeT>::node_iterator&
PeriodicSCLattice<nodeT>::node_iterator::to_neighbour(int i)
{
  switch(i) {
  case 0: return to(&N());
  case 1: return to(&E());
  case 2: return to(&S());
  case 3: return to(&W());
  case 4: return to(&U());
  case 5: return to(&D());
  }
}


template <typename nodeT,typename Function>
struct implement_for_each_neighbour<nodeT,Function,typename PeriodicSCLattice<nodeT>::node_iterator> {
  inline static Function fen(typename PeriodicSCLattice<nodeT>::node_iterator n,Function f)
  {
    f(n.N());
    f(n.E());
    f(n.S());
    f(n.W());
    f(n.U());
    f(n.D());
    return f;
  }
} ;

template <typename nodeT,typename Function>
struct implement_for_each_bond<nodeT,Function,PeriodicSCLattice<nodeT>> {
  inline static Function fen(PeriodicSCLattice<nodeT> &lat,Function f)
  {
    for (auto iter=lat.begin(); iter!=lat.end(); ++iter) {
      f(*iter,iter.W());
      f(*iter,iter.S());
      f(*iter,iter.D());
    }
    return f;
  }
} ;

} /* namespace */

#endif /* LATTICE3D_HH */
