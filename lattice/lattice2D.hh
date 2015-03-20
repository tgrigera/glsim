/*
 * lattice2D.hh -- Lattices in two dimensions
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

#ifndef _LATTICE2D_HH_
#define _LATTICE2D_HH_

#include <graph.hh>

/* @T
\section{Periodic square lattice}

@c */

template <typename nodeT>
class PeriodicSQLattice;

template <typename nodeT>
class PeriodicSQLattice : public GraphBase<nodeT> {
public:
  PeriodicSQLattice(int Lx,int Ly);
  PeriodicSQLattice(const PeriodicSQLattice &);
  ~PeriodicSQLattice();
  PeriodicSQLattice& operator=(const PeriodicSQLattice&);

  // Node access

  int size_x() const {return Lx;}
  int size_y() const {return Ly;}
  nodeT&       operator()(int x,int y) {return nodes[lpos(x,y)];}
  const nodeT& operator()(int x,int y) const {return nodes[lpos(x,y)];}
  ptrdiff_t    id(int x,int y) const {return lpos(x,y);}
  ptrdiff_t    id_protected(int x,int y) const {
    return lpos(mymod(x,Lx),mymod(y,Ly));}
  using GraphBase<nodeT>::id;

  // I/O

  void write(std::ofstream&);
  void read(std::ifstream&);

  // Public utility
  int pdiff(int i,int j,int L) const;
  int pdiffx(int i, int j) const {return pdiff(i,j,Lx);}
  int pdiffy(int i, int j) const {return pdiff(i,j,Ly);}


  // Iterators

  class neighbour_iterator;
  class node_iterator;
  node_iterator  node_begin()
  {return node_iterator(*this);}
  nodeT*         node_end()
  {return GraphBase<nodeT>::node_end();}
  neighbour_iterator neighbour_begin(ptrdiff_t id)
  {return neighbour_begin(this->begin()+id);}
  neighbour_iterator neighbour_begin(nodeT* n)
  {return neighbour_iterator(*this,n);}
  neighbour_iterator neighbour_begin(ptrdiff_t i,ptrdiff_t j)
  {return neighbour_iterator(*this,i,j);}
  int neighbour_end(ptrdiff_t id)
  {return 4;}
  int neighbour_end(nodeT* n)
  {return 4;}

private:
public:
  int  Lx,Ly;
  int  *Ndis,*Sdis,*Edis,*Wdis;

  // This transforms coordinates to a valid node id
  ptrdiff_t lpos(ptrdiff_t i,ptrdiff_t j) const {return Ly*i + j;}
  void posl(ptrdiff_t l,ptrdiff_t &i,ptrdiff_t &j);

  using GraphBase<nodeT>::nodes;

  void load_neighbour_distance_vectors();

  int mymod(int a,int b) {int r=a%b; return (r>=0) ? r : r+b;}

} ;

/* @T

\paragraph{Construction}

@c */

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

template<typename nodeT>
PeriodicSQLattice<nodeT>::~PeriodicSQLattice()
{
  delete[] Ndis;
  delete[] Sdis;
  delete[] Wdis;
  delete[] Edis;
}

/* @T
\section{whatever}
@c */


template<typename nodeT>
inline void PeriodicSQLattice<nodeT>::posl(ptrdiff_t l,ptrdiff_t &i,ptrdiff_t &j)
{
  i=l/Ly;
  j=l-i*Ly;
}

template<typename nodeT>
void PeriodicSQLattice<nodeT>::write(std::ofstream& f)
{
  GraphBase<nodeT>::write(f);
  f.write((char*) &Lx,sizeof(Lx));
  f.write((char*) &Ly,sizeof(Ly));
}

template<typename nodeT>
void PeriodicSQLattice<nodeT>::read(std::ifstream& f)
{
  GraphBase<nodeT>::read(f);
  f.read((char*) &Lx,sizeof(Lx));
  f.read((char*) &Ly,sizeof(Ly));
  load_neighbour_distance_vectors();
}

template<typename nodeT>
inline int PeriodicSQLattice<nodeT>::pdiff(int i,int j,int L) const
{
  return i-j-L*rint((double)(i-j)/L);
}

/* @T
\section{Iterators}
@c*/

template <typename nodeT>
class PeriodicSQLattice<nodeT>::neighbour_iterator :
  public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:

  typedef PeriodicSQLattice<nodeT> Graph_t;

  neighbour_iterator(Graph_t& lat,int i,int j) :
    nn(nlist)
  {
    nodeT* n=Graph_t::begin()+lat.lpos(i,j);
    load_nlist(lat,n,i,j);
  }

  neighbour_iterator(Graph_t& lat,nodeT* n) :
    nn(nlist)
  {
    ptrdiff_t i,j;
    lat.posl(lat.id(n),i,j);
    load_nlist(lat,n,i,j);
  }

  neighbour_iterator(const neighbour_iterator &ni)
  {memcpy(nlist,ni.nlist,4*sizeof(nodeT*));
    nn=nlist+(ni.nn-ni.nlist);}

  neighbour_iterator& operator=(const neighbour_iterator& ni)
  {memcpy(nlist,ni.nlist,4*sizeof(nodeT*));
    nn=nlist+(ni.nn-ni.nlist);}

  bool operator==(const typename GraphBase<nodeT>::neighbour_iterator &i)
  {return nn==i.nn;}

  bool operator!=(int p)
  {return nn!=nlist+p;}

  typedef typename GraphBase<nodeT>::iterator iterator;

  operator iterator() {return *nn;}

  nodeT& operator*() const {return **nn;}
  nodeT* operator->() const {return *nn;}

  neighbour_iterator& operator++() {++nn; return *this;}
  neighbour_iterator& operator++(int) {neighbour_iterator c=*this; ++nn; return *this;}
  neighbour_iterator& operator--() {nn--; return *this;}
  neighbour_iterator& operator--(int) {neighbour_iterator c=*this; --nn; return *this;}

private:
  nodeT **nn;
  nodeT *nlist[4];

  void load_nlist(Graph_t& lat,nodeT* n,int i,int j)
  {
    nlist[0]=n+lat.Ndis[j];
    nlist[1]=n+lat.Sdis[j];
    nlist[2]=n+lat.Edis[i];
    nlist[3]=n+lat.Wdis[i];
  }
} ;



template <typename nodeT>
class PeriodicSQLattice<nodeT>::node_iterator :
  public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:
  typedef PeriodicSQLattice<nodeT> Graph_t;

  node_iterator(PeriodicSQLattice<nodeT> &l,ptrdiff_t i_=0,ptrdiff_t j_=0) :
    lat(l), i(i_), j(j_), n(l.nodes+l.id(i,j)) {}

  node_iterator(const node_iterator &i) :
  lat(i.lat), i(i.i), j(i.j), n(i.n) {}

  node_iterator& operator=(const node_iterator &it)
  {if (this==&it) return *this; lat=it.lat; i=it.i; j=it.j; n=it.n; return *this;}

  bool operator==(const node_iterator &i)
  {return n==i.n;}

  bool operator!=(const node_iterator &i)
  {return n!=i.n;}

  bool operator==(nodeT* p)
  {return n==p;}

  bool operator!=(nodeT* p)
  {return n!=p;}

  nodeT* N() const {return n+lat.Ndis[j];}
  nodeT* S() const {return n+lat.Sdis[j];}
  nodeT* E() const {return n+lat.Edis[i];}
  nodeT* W() const {return n+lat.Wdis[i];}

  nodeT& operator*() const {return *n;}
  nodeT* operator->() const {return n;}

  node_iterator& operator++() {++n; i+=j/(lat.Ly-1); j+=lat.Ndis[j]; return *this;}
  node_iterator operator++(int) {node_iterator c=*this; ++(*this); return c;}

  node_iterator& operator--() {--n; j+=lat.Sdis[j]; i-=j/(lat.Nx-1);   return *this;}
  node_iterator operator--(int) {node_iterator c=*this; --(*this); return c;}

  ptrdiff_t  ii() {return i;}
  ptrdiff_t  jj() {return j;}
  
  // Above methods are required by a standard iterator; we add the following

  // node_iterator& to(int n) {nn=graph.nodes+nn; return *this;}
  // node_iterator& to_neighbour(int n) {nn=graph.node_neighbours(nn)+n; return *this;}
  // nodeT* neighbour(int i) const {return graph.node_neighbours(nn)+i;}

// protected:
public:
  Graph_t&  lat;
  nodeT*    n;
  ptrdiff_t i,j;
} ;


template <typename nodeT,typename Function>
inline Function for_each_neighbour(typename PeriodicSQLattice<nodeT>::node_iterator& n,Function f)
{
  f(n.N());
  f(n.S());
  f(n.E());
  f(n.W());
}

#endif /* _LATTICE2D_HH_ */
