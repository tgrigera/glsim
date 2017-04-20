/*
 * Poisson.hh -- Poisson graphs
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

#ifndef POISSON_HH
#define POISSON_HH

#include <limits>

#include "random.hh"
#include "olconfiguration.hh"
#include "nneighbours.hh"
#include "graph.hh"

namespace glsim {

/******************************************************************************
 *
 * Poisson lattice in 3-d (metric, periodic)
 *
 */

/** \ingroup lattice

    \brief Metric Poisson graph in 3-d

    This class implements a metric Poisson graph in 3-d, i.e. a graph
    formed by taking N points randomly and independently placed within
    a parallelepiped and connecting those lying within an Euclidean
    distance less than a given cut-off.
*/
template <typename nodeT>
class MetricPoisson3D : public GraphBase<nodeT> {
public:
  /// \name Construction, copy and destruction
  ///@{
  MetricPoisson3D(int N,double box[],double irange,bool periodic=true);
  MetricPoisson3D(const MetricPoisson3D &);
  MetricPoisson3D& operator=(const MetricPoisson3D&);

  ///@}
  /// \name Lattice information
  ///@{
  int size() const {return conf.N;}
  int neighbour_size(id_t id) const {return nlist.Nneighbours(id);}

  ///@}
  /// \name Node access
  ///@{
  using GraphBase<nodeT>::id;
  using GraphBase<nodeT>::data;
  
  ///@}
  ///\name Disk i/o
  ///@{
  void read(std::istream&);  ///< Read graph from stream
  void write(std::ostream&); ///< Write graph to stream

  ///@}

  ///@}
  ///\name Utility
  ///@{
  const OLconfiguration& oconf() const {return conf;}  ///< Const reference to OLconfiguration object for node positions

  ///@}
  /// \name Iterators
  ///@{
  class   node_iterator;

  node_iterator      begin() {return node_iterator(*this);}
  nodeT*             end()   {return GraphBase<nodeT>::end();}
  ///@}

private:
  double                 intrange;
  OLconfiguration        conf;
  NeighbourList_subcells nlist;

  void init(int,double[]);

  template <typename nodeTa,typename Function,typename Iterator>
  friend struct implement_for_each_neighbour;
  template <typename nodeTa,typename Function,typename Graph>
  friend struct implement_for_each_bond;
} ;

/*******************************************************************************
 * Construction, copy and destruction
 *
 */

/**
   Create and build graph.

  \param N      Number of nodes
  \param box    Size of parallelepiped containing the nodes
  \param irange Interaction (cut-off) range
  \param periodic  Whether to compute the distances using periodic boundary conditions
 
 */
template <typename nodeT> inline
MetricPoisson3D<nodeT>::MetricPoisson3D(int N,double box[],double irange,
					bool periodic) :
  GraphBase<nodeT>(N,-1,false),
  intrange(irange),
  nlist(irange,0.)
{
  if (!periodic)
    box[0]=box[1]=box[2]=std::numeric_limits<double>::max();
  init(N,box);
}

template <typename nodeT> inline
MetricPoisson3D<nodeT>::MetricPoisson3D(const MetricPoisson3D& g) :
  GraphBase<nodeT>(g),
  intrange(g.intrange),
  conf(g.conf),
  nlist(g.intrange,0.)
{
  nlist.rebuild(conf,intrange);
}  

template<typename nodeT>
MetricPoisson3D<nodeT>& MetricPoisson3D<nodeT>::
operator=(const MetricPoisson3D<nodeT>& g)
{
  if (this==&g) return *this;
  GraphBase<nodeT>::operator=(g);
  intrange=g.intrange;
  conf=g.conf;
  nlist.rebuild(conf,intrange);
}

template <typename nodeT>
void MetricPoisson3D<nodeT>::init(int N,double boxl[])
{
  conf.N=N;
  conf.time=0;
  conf.step=0;
  conf.box_angles[0]=conf.box_angles[1]=conf.box_angles[2]=90.;
  conf.box_length[0]=boxl[0];
  conf.box_length[1]=boxl[1];
  conf.box_length[2]=boxl[2];
  conf.r=new double[conf.N][3];
  conf.v=conf.a=0;
  conf.id=conf.type=conf.flags=0;

  Uniform_real ran(0,1);
  for (int i=0; i<conf.N; ++i) {
    conf.r[i][0]=conf.box_length[0]*ran();
    conf.r[i][1]=conf.box_length[1]*ran();
    conf.r[i][2]=conf.box_length[2]*ran();
  }
  nlist.rebuild(conf,intrange);
}

/******************************************************************************
 *
 * Utilitty and i/o
 *
 */
template<typename nodeT>
void MetricPoisson3D<nodeT>::write(std::ostream& f)
{
  GraphBase<nodeT>::write(f);
  f.write((char*) &intrange,sizeof(intrange));
  f.write((char*) &conf.N,sizeof(conf.N));
  f.write((char*) &conf.time,sizeof(conf.time));
  f.write((char*) &conf.step,sizeof(conf.step));
  f.write((char*) &conf.box_length,sizeof(conf.box_length));
  f.write((char*) &conf.box_angles,sizeof(conf.box_angles));
  f.write((char*) conf.r,3*conf.N*sizeof(double));
}

template<typename nodeT>
void MetricPoisson3D<nodeT>::read(std::istream& f)
{
  GraphBase<nodeT>::read(f);
  f.read((char*) &intrange,sizeof(intrange));
  f.read((char*) &conf.N,sizeof(conf.N));
  f.read((char*) &conf.time,sizeof(conf.time));
  f.read((char*) &conf.step,sizeof(conf.step));
  f.read((char*) &conf.box_length,sizeof(conf.box_length));
  f.read((char*) &conf.box_angles,sizeof(conf.box_angles));
  delete[] conf.r;
  conf.r=new double[conf.N][3];
  f.read((char*) conf.r,3*conf.N*sizeof(double));
  nlist.rebuild(conf,intrange);
}

/******************************************************************************
 *
 * Iterator
 *
 */
/** \brief A bidirectional iterator for the 3-d metric Poisson graph

 */
template <typename nodeT>
class MetricPoisson3D<nodeT>::node_iterator :
  public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:
  typedef MetricPoisson3D<nodeT> Graph_t;  ///< The type of the graph we belong to

  ///\name Construction, copy and comparison
  //@{
  node_iterator(MetricPoisson3D<nodeT> &l,nodeT* n);
  ///< Create and move to n
  node_iterator(MetricPoisson3D<nodeT> &l,id_t i=0);
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
  {++n; return *this;}

  node_iterator operator++(int) {node_iterator c=*this; ++(*this); return c;}

  node_iterator& operator--()
  {--n; return *this;}

  node_iterator operator--(int) {node_iterator c=*this; --(*this); return c;}
  //@}

  ///\name Extra methods (including neighbour access)
  ///@{
  node_iterator& to(nodeT* n);  ///< Move to node given by pointer.
  node_iterator& to(id_t i);  ///< Move to node given by id.

  int             neighbour_size() const {return grph.nlist.Nneighbours(grph.id(n));}
  ///< Coordination number of node (always 6).
  node_iterator& to_neighbour(int i); ///<Move to `i`th neighbour (numbers as in neighbour()).

  nodeT& neighbour(int i) const;
  ///< Access neighbour by number

  operator nodeT*() const {return n;}

  ///@}

private:
  Graph_t&  grph;
  nodeT*    n;

  template <typename nodeTa,typename Function,typename Iterator>
  friend struct implement_for_each_neighbour;
} ;

/*******************************************************************************
 *
 * Construction, destruction and comparison
 *
 */

template <typename nodeT> inline MetricPoisson3D<nodeT>::node_iterator::
node_iterator(MetricPoisson3D<nodeT> &g ,nodeT* n_) :
  grph(g),
  n(n_)
{}

template <typename nodeT> inline MetricPoisson3D<nodeT>::node_iterator::
node_iterator(MetricPoisson3D<nodeT> &g ,id_t i) :
  grph(g),
  n(g.data()+i)
{}

template <typename nodeT> inline MetricPoisson3D<nodeT>::node_iterator::
node_iterator(const node_iterator &i) :
 grph(i.grph),
 n(i.n)
{}

template <typename nodeT> inline
typename MetricPoisson3D<nodeT>::node_iterator&
MetricPoisson3D<nodeT>::node_iterator::operator=(const node_iterator &it)
{
  if (this==&it) return *this;
  grph=it.grph;
  n=it.n;
  return *this;
}

/*
 * Movement and neighbours
 *
 */
template <typename nodeT> inline
typename MetricPoisson3D<nodeT>::node_iterator&
MetricPoisson3D<nodeT>::node_iterator::to(nodeT* node)
{
  n=node;
  return *this;
}

template <typename nodeT> inline
typename MetricPoisson3D<nodeT>::node_iterator&
MetricPoisson3D<nodeT>::node_iterator::to(id_t i)
{
  return to(grph.data()+i);
}
 
template <typename nodeT> inline
nodeT& MetricPoisson3D<nodeT>::node_iterator::neighbour(int i) const
{
  return *(grph.data()+grph.nlist.neighbour_list(grph.id(n)).at(i));
}

template <typename nodeT> inline
typename MetricPoisson3D<nodeT>::node_iterator&
MetricPoisson3D<nodeT>::node_iterator::to_neighbour(int i)
{
  n=grph.data()+grph.nlist.neighbour_list(grph.id(n)).at(i);
}

template <typename nodeT,typename Function>
struct implement_for_each_neighbour<nodeT,Function,typename MetricPoisson3D<nodeT>::node_iterator> {
  inline static Function fen(typename MetricPoisson3D<nodeT>::node_iterator n,Function f)
  {
    for (auto pn=n.grph.nlist.neighbours_begin(n.grph.id(n.n)), end=n.grph.nlist.neighbours_end(n.grph.id(n.n));
	 pn!=end; ++pn)
      f(n.grph[*pn]);
    return f;
  }
} ;

template <typename nodeT,typename Function>
struct implement_for_each_bond<nodeT,Function,MetricPoisson3D<nodeT>> {
  inline static Function fen(MetricPoisson3D<nodeT> &lat,Function f)
  {
    for (auto p = lat.nlist.pairs_begin(), end=lat.nlist.pairs_end(); p!=end; ++p) {
#ifdef DEBUG
      if (lat.conf.distancesq(p->first,p->second)>lat.intrange)
	throw Runtime_error("This should never happen",HERE);
#endif
      f(lat[p->first],lat[p->second]);
    }
    return f;
  }
} ;

}  /* glsim */

#endif /* POISSON_HH */
