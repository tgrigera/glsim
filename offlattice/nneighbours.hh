/*
 * nneighbours.hh -- structures to find nearest neighbours
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

#include <utility>
#include <vector>

#include "olconfiguration.hh"

#ifndef NNEIGHBOURS_HH
#define NNEIGHBOURS_HH

namespace glsim {

/*****************************************************************************
 *
 * Metric nearest neighbours
 *
 */

/** \class MetricNearestNeighbours
    \ingroup OfflatticeINT

NOTE!!!! Algorithms for nearest neighbours are to be thoght as
providing __candidates__for NN.  For most algorithms here, pairs
returned are not guaranteed to be within cutoff, the algorithms here
mereley reduce the candidates list.  This is a feature, because in
this way the structure need not be rebuilt at every step (which may be
costly).  Candidates are guaranteed to appear only once.  Check
Interactions for use.

MetricNearestNeighbours a class that illustrates the interface of the
classes that find the metric nearest neighbours (i.e. neighbours
defined as those lying at a distance less than some value, contrast
with TopologicalNearestNeighbours).  It is just an illustration
because the different implementations of metric nearest-neighbour
search do not derive from this.  Instead, they are inteded to be used
in a template class like interactions, with a STL-container-like
interface.  In particular the functions pairs_begin() etc return an
appropriate iterator, but its exact type is not known here, so the
implementations don't inherit from this.  Type erasure
(http://www.artima.com/cppsource/type_erasure.html) could help here,
but I'm not yet convinced it's worth the effort.

This class simply returns all pairs as candidates, serving only as
interface illustration and for benchmarking.

*/
class MetricNearestNeighbours {
public:
  typedef std::pair<int,int> Pair;
  
  MetricNearestNeighbours(double rc_) : rc(rc_) {rcsq=rc*rc;}
  virtual ~MetricNearestNeighbours() {}

  /// Build lists from scratch for given conf
  virtual void rebuild(OLconfiguration& c,double rc_=-1)
  {conf=&c; if (rc_>0) rc=rc_;}
  /// Inform of change in configuration, will try to update lists
  /// assuming particles have not moved much, may rebuild everything
  /// from scratch
  virtual void update(double a=0) {}

  double  cutoff()   const {return rc;}
  double  cutoffsq() const {return rcsq;}

  class NeighbourIterator;
  NeighbourIterator neighbours_begin(int i);
  int               neighbours_end(int i) {return -1;}

  class PairIterator;
  PairIterator      pairs_begin();
  int               pairs_end() {return -1;}

  template <typename Function,typename NeighboursT>
  friend struct implement_for_each_pair;

private:
  OLconfiguration* conf;
  double           rc,rcsq;
} ;

class MetricNearestNeighbours::NeighbourIterator :
    public std::iterator<std::input_iterator_tag,int>
{
public:
  NeighbourIterator(MetricNearestNeighbours& NN_,int particle_) :
    NN(NN_),
    particle(particle_)
  { n = particle==0 ? 1 : 0; }

  bool operator==(int i) const
  {return n>=NN.conf->N;}

  bool operator!=(int) const
  {return n<NN.conf->N;}

  const int& operator*() const {return n;}
  int const * operator->() const {return &n;}

  NeighbourIterator& operator++()
  { while (++n==particle); }

  NeighbourIterator operator++(int)
  {NeighbourIterator c=*this; ++(*this); return c;}

private:
  MetricNearestNeighbours& NN;
  int  particle;   ///< Which particle are we tracking neighbours to
  int  n;          ///< Current neighbour
} ;

inline MetricNearestNeighbours::NeighbourIterator
MetricNearestNeighbours::neighbours_begin(int i)
{
  return NeighbourIterator(*this,i);
}

class MetricNearestNeighbours::PairIterator :
    public std::iterator<std::input_iterator_tag,int>
{
public:
  PairIterator(MetricNearestNeighbours &NN_) :
    NN(NN_), pair(0,1) {}

  bool operator==(const PairIterator& ni) const
  {return pair==ni.pair;}
  bool operator!=(const PairIterator& ni) const
  {return ! operator==(ni);}

  bool operator==(int)  const
  {return pair.first>=NN.conf->N-1;}
  bool operator!=(int)  const
  {return pair.first<NN.conf->N-1;}

  const Pair& operator*() const {return pair;}
  Pair const * operator->() const {return &pair;}

  PairIterator& operator++()
  {++pair.second; if (pair.second==NN.conf->N) {++pair.first; pair.second=pair.first+1;}}

  PairIterator operator++(int)
  {PairIterator c=*this; ++(*this); return c;}

private:
  MetricNearestNeighbours NN;
  Pair                    pair;
} ;

inline MetricNearestNeighbours::PairIterator
MetricNearestNeighbours::pairs_begin()
{
  return PairIterator(*this);
}

/** \class NeighbourList_naive
    \ingroup OfflatticeINT

    This provides a naive implementation of a neighbour list.
    Although the list is built by iterating through all possible
    pairs, if the simulation step provides information on the maximum
    displacement, this brings noticeable speed improvements at more
    than 200 particles, at least for potentials with relatively short
    cut-offs (such as the repulsive Lennard-Jones).
*/
class NeighbourList_naive {
public:
  typedef std::pair<int,int> Pair;

  NeighbourList_naive(double rc,double delta_r=-1);
  void rebuild(OLconfiguration&,double rc=-1);
  void update(double maxdisp);

  double  cutoff() const {return rc;}
  double  cutoffsq() const {return rcsq;}

  std::vector<Pair>::iterator pairs_begin() {return pairs.begin();}
  std::vector<Pair>::iterator pairs_end() {return pairs.end();}
  std::vector<int>::iterator neighbours_begin(int i) {return neighbours[i].begin();}
  std::vector<int>::iterator neighbours_end(int i) {return neighbours[i].end();}

  template <typename Function,typename NeighboursT>
  friend struct implement_for_each_pair;

private:
  OLconfiguration*  conf;
  std::vector<Pair> pairs;
  std::vector<std::vector<int>>  neighbours;

  double rc,rcsq,rdsq,delta_r;
  double accum_maxdisp;
} ;


/** \class Subcells
    \ingroup OfflatticeINT

    This class implements neighbour lookup throug subcells [Allen].

*/
class Subcells {
public:
  typedef std::pair<int,int> Pair;

  Subcells(double rc,double delta_r=0);
  /// Build lists from scratch for given configuration. Will check for
  /// changes in box size and number of particles.
  void rebuild(OLconfiguration&,double rc=-1);
  /// Inform of change in configuration, will try to update lists
  /// assuming particles have not moved much, may rebuild everything
  /// from scratch.  Assumes number of particles and box size has not
  /// changed.
  void update(double maxdisp);
  ~Subcells()       {clear_lists();};

  double cutoff() const {return rc;}
  double cutoffsq() const {return rcsq;}

  class NeighbourIterator;
  NeighbourIterator neighbours_begin(int i);
  NeighbourIterator neighbours_end(int i);

  class PairIterator;
  PairIterator      pairs_begin();
  PairIterator      pairs_end();

  void  init_cells(int,double[]);
  void  clear_lists();
  int   num_subcells() const;
  int   first_particle(int cell) const;
  int   next_particle(int particle) const;
  int   neighbour(int cell,int n) const;
  int   which_subcell(int particle) const;

  template <typename Function,typename NeighboursT>
  friend struct implement_for_each_pair;

private:
  OLconfiguration* conf;
  double rc,rcsq,delta_r;
  double accum_maxdisp;

  int    m[3];                   ///< Number of subcells in each axis
  int    nscell;                 ///< Total number of subcells in box
  int    nparticles;             ///< Number of particles
  double scell[3];               ///< Size of subcells in each axis 
  int    *subcell;               ///< Head particle of each subcel
  int    *llist;                 ///< Linked list
  int    *WhichSubcell;          ///< To find subcell from particle index
  int    *subcelln;              ///< Subcell neighbour list
  double  boxl[3];               ///< Box size used to build subcell structure

  int   icell(int ix,int iy,int iz);
  int   eicell(int ix,int iy,int iz);
} ;


/**
 * \param delta_r  Cell size will be rc + delta_r
 *
 * delta_r Should be equal to DR in MC and 0 in for MD
 *
 */
inline Subcells::Subcells(double rc_,double delta_r_) :
  rc(rc_), delta_r(delta_r_),
  accum_maxdisp(0),
  nscell(0),nparticles(0),
  subcell(0),llist(0),
  WhichSubcell(0),
  subcelln(0)
{
  m[0]=m[1]=m[2]=0;
}

inline int Subcells::icell(int ix,int iy,int iz)
{
  return ix*m[1]*m[2] + iy*m[2] + iz;
}

inline int Subcells::eicell(int ix,int iy,int iz)
{
  return ((ix+m[0])%m[0])*m[1]*m[2] + ((iy+m[1])%m[1])*m[2] +
    (iz+m[2])%m[2];
}

inline int Subcells::num_subcells() const
{
  return nscell;
}

inline int Subcells::first_particle(int cell) const
{
  return subcell[cell];
}

inline int Subcells::next_particle(int particle) const
{
  return llist[particle];
}

inline int Subcells::neighbour(int cell,int n) const
{
  return subcelln[cell*27+n];
}

inline int Subcells::which_subcell(int particle) const
{
  return WhichSubcell[particle];
}

//
// The iterator for neighbours
//
class Subcells::NeighbourIterator :
    public std::iterator<std::input_iterator_tag,int>
{
public:
  NeighbourIterator(Subcells &sub,int particle,bool end=false);
  NeighbourIterator(const NeighbourIterator& ni);

  bool operator==(const NeighbourIterator& ni) const;
  bool operator==(int) const {return n==-1;}
  bool operator!=(const NeighbourIterator& ni) const
  {return ! operator==(ni);}
  bool operator!=(int) const {return n!=-1;}

  const int& operator*() const {return n;}
  int const * operator->() const {return &n;}

  NeighbourIterator& operator++();

  NeighbourIterator operator++(int)
  {NeighbourIterator c=*this; ++(*this); return c;}

private:
  Subcells& scell; ///< Our subcell structure
  int  particle;   ///< Which particle are we tracking neighbours to
  int  icell;      ///< Cell index for the cell of our particle
  int  n;          ///< Current neighbour
  int  ncell;      ///< Current neighbouring cell (-1 is self, max 26)
} ;

inline Subcells::NeighbourIterator::NeighbourIterator(Subcells& sub,int part_,bool end) :
  scell(sub),
  particle(part_),
  ncell(-1)
{
  if (end) {n=-1; return;}
  icell=scell.which_subcell(particle);;
  n=scell.first_particle(icell);
  if (n==particle) ++(*this);
}

inline Subcells::NeighbourIterator::
NeighbourIterator(const Subcells::NeighbourIterator &ni) :
  scell(ni.scell),
  particle(ni.particle),
  icell(ni.icell),
  n(ni.n),
  ncell(ni.ncell)
{}
  
inline bool Subcells::NeighbourIterator::
operator==(const Subcells::NeighbourIterator& ni) const
{
  return particle==ni.particle && n==ni.n;
}

class System_too_small : public Runtime_error {
public:
  explicit System_too_small(const std::string& dir,double actual,double min,
			    const Source_context &c=Source_context()) :
    Runtime_error("Linear size in " + dir + " direction too small ("
		  + std::to_string(actual) + ", minimum "+std::to_string(min) +")",
		  c)
  {}
} ;

//
// The iterator for pairs
//
class Subcells::PairIterator :
    public std::iterator<std::input_iterator_tag,int>
{
public:
  PairIterator(Subcells &sub,bool end=false);
  PairIterator(const PairIterator& ni);

  bool operator==(const PairIterator& ni) const;
  bool operator!=(const PairIterator& ni) const
  {return ! operator==(ni);}

  const Pair& operator*() const {return pair;}
  Pair const * operator->() const {return &pair;}

  PairIterator& operator++();

  PairIterator operator++(int)
  {PairIterator c=*this; ++(*this); return c;}

private:
  Subcells& scell; ///< Our subcell structure
  Pair      pair;  ///< Current pair of neighbours
  int       icell; ///< Cell index of first member of pair
  int       ncell; ///< Neighbour (cell) number
} ;

inline Subcells::PairIterator::PairIterator(Subcells& sub,bool end) :
  scell(sub),
  pair(-1,-1),
  icell(-1),
  ncell(-1)
{
  if (end) return;
  while (pair.first<0)
    pair.first=scell.first_particle(++icell);
  pair.second=pair.first;
  ++(*this);  // let operator++ find the first valid jparticle (which
	      // may be in a different cell
}

inline Subcells::PairIterator::
PairIterator(const Subcells::PairIterator &ni) :
  scell(ni.scell),
  pair(ni.pair),
  icell(ni.icell),
  ncell(ni.ncell)
{}
  
inline bool Subcells::PairIterator::
operator==(const Subcells::PairIterator& ni) const
{
  return &scell==&(ni.scell) && pair==ni.pair;
}

//
// The functions returning iterators
//

inline Subcells::NeighbourIterator Subcells::neighbours_begin(int i)
{
  return NeighbourIterator(*this,i);
}

inline Subcells::NeighbourIterator Subcells::neighbours_end(int i)
{
  return NeighbourIterator(*this,i,true);
}

inline Subcells::PairIterator Subcells::pairs_begin()
{
  return PairIterator(*this);
}

inline Subcells::PairIterator Subcells::pairs_end()
{
  return PairIterator(*this,true);
}

/** \class NeighbourList_subcells
    \ingroup OfflatticeINT

    This class provides metric nearest neighbour candidates through a
    pair list.  The pair list is built using subcells (class
    Subcells).  The provided candidates will be separated by a
    distance at most rc + delta_r.
*/
class NeighbourList_subcells {
public:
  typedef std::pair<int,int> Pair;

  NeighbourList_subcells(double rc,double delta_r=-1);
  void rebuild(OLconfiguration&,double rc=-1);
  void update(double maxdisp);

  double cutoff() const {return rc;}
  double cutoffsq() const {return rcsq;}

  std::vector<Pair>::iterator pairs_begin() {return pairs.begin();}
  std::vector<Pair>::iterator pairs_end() {return pairs.end();}
  std::vector<int>::iterator neighbours_begin(int i) {return neighbours[i].begin();}
  std::vector<int>::iterator neighbours_end(int i) {return neighbours[i].end();}

  template <typename Function,typename NeighboursT>
  friend struct implement_for_each_pair;

private:
  OLconfiguration   *conf;
  std::vector<Pair> pairs;
  std::vector<std::vector<int>>  neighbours;

  double    rc,rcsq,rdsq,delta_r;
  double    accum_maxdisp;

  Subcells  SUBC;
} ;

/*****************************************************************************
 *
 * Functions to loop through all nearest neigbhours, with the static
 * member function trick for partial specialization
 *
 */
template <typename Function,typename NeighboursT>
struct implement_for_each_pair {
  static Function for_each_pair(NeighboursT& NN,Function f)
  {
    for (auto p = NN.pairs_begin(), end=NN.pairs_end(); p!=end; ++p) {
      double dsq=NN.conf->distancesq(p->first,p->second);
      if (dsq<=NN.cutoffsq()) f(p->first,p->second,dsq);
    }
    return f;
  }
} ;

template <typename Function>
struct implement_for_each_pair<Function,glsim::MetricNearestNeighbours> {
  static Function for_each_pair(glsim::MetricNearestNeighbours& NN,Function f)
  {
    for (int i=0; i<NN.conf->N-1; ++i)
      for (int j=i+1; j<NN.conf->N; ++j) {
	double dsq=NN.conf->distancesq(i,j);
	if (dsq<=NN.cutoffsq()) f(i,j,dsq);
      }
    return f;
  }
} ;

template <typename Function>
struct implement_for_each_pair<Function,glsim::Subcells> {
  static Function for_each_pair(glsim::Subcells& NN,Function f)
  {
    for (int isc=0; isc<NN.num_subcells(); isc++) /* Loop over subcells */ 
      for (int n=NN.first_particle(isc); n>=0; n=NN.next_particle(n)) { /* and particles */                                                                       
	/* find pairs in the same subcell */ 
	for (int m=NN.next_particle(n); m>=0; m=NN.next_particle(m)) { 
	  double dsq=NN.conf->distancesq(n,m); 
	  if (dsq<NN.cutoffsq()) f(n,m,dsq);
	} 
	/* find pairs in neighbouring subcells */ 
	for (int nn=0; nn<13; nn++)
	  for (int m=NN.first_particle(NN.neighbour(isc,nn)); m>=0; m=NN.next_particle(m) ) {                                                                    
	    double dsq=NN.conf->distancesq(n,m); 
	    if (dsq<NN.cutoffsq()) f(n,m,dsq); 
	  } 
      } 
    return f;
  }
} ;

template <typename Function,typename NeighboursT>
Function for_each_pair(NeighboursT& NN,Function f)
{
  implement_for_each_pair<Function,NeighboursT>::for_each_pair(NN,f);
}

/*****************************************************************************
 *
 * Topological nearest neighbours
 *
 */

/** \class TopologicalNearestNeighbours
    \ingroup OfflatticeINT

    TopologicalNearestneighbours is an abstract base class that serves
as front end to different implementations of topological
nearest-neighbour search (i.e. neighbourhoods with a fixed number of
nearest neighbours, contrast with MetricNearestNeighbours), thus
leading to potentially non-symmetric interaction matrices.  It
supports rebuilding on demand or updates (which may imply complete
rebuilding) controlled by the maximum modulus of particle
displacements.
*/
class TopologicalNearestNeighbours {
public:
  TopologicalNearestNeighbours(int Nnearest_) : Nnearest(Nnearest_) {}
  virtual ~TopologicalNearestNeighbours() {}

  /// Build list from scratch (depending on the actual algorithm, this
  /// may not actually build a list, but will reinitialize everything
  /// for a new configuration).
  virtual void rebuild(OLconfiguration&,int n)=0;
  /// Inform of change in configuration, will try to update lists
  /// assuming particles have not moved much, may rebuild everything
  /// from scratch
  virtual void update(double)=0;

  virtual std::vector<int>::iterator neighbours_begin(int i)=0;
  virtual std::vector<int>::iterator neighbours_end(int i)=0;

protected:
  int    Nnearest;
} ;

/** \class TopologicalNeighbours_naive
    \ingroup OfflatticeINT

    This produces pairs of nearests neighbours, found with the usual
    Euclidean distance, but for each particle finds exactly N nearest
    neighbours,  Implementation is naive (linear search of neighbours
    for each particle).  It is rather slow and meant mainly as a check
    for implementations with better algorithms.
*/
class TopologicalNeighbours_naive : public TopologicalNearestNeighbours {
public:
  TopologicalNeighbours_naive(int Nnearest) :
    TopologicalNearestNeighbours(Nnearest) {};
  void rebuild(OLconfiguration&,int Nnearest=-1);
  // maybe this makes no sense:
  void update(double maxdisp) {rebuild(*conf);}

  std::vector<int>::iterator neighbours_begin(int i) {return neighbours[i].begin();}
  std::vector<int>::iterator neighbours_end(int i) {return neighbours[i].end();}

private:
  OLconfiguration* conf;
  struct ndist {int n; double dsq; ndist(int n_,double d_) : n(n_), dsq(d_){} };
  std::vector<std::vector<int>>  neighbours;
} ;


} /* namespace */

#endif /* NNEIGHBOURS_HH */
