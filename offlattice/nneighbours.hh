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

#include <vector>

#include "olconfiguration.hh"

#ifndef NNEIGHBOURS_HH
#define NNEIGHBOURS_HH

namespace glsim {

/** \class NearestNeighbours
    \ingroup OfflatticeINT

    NearestNeighbours is an abstract base class that serves as front
end to different implementations of nearest-neighbour search.  It
supports rebuilding on demand or updates (which may imply complete
rebuilding) controlled by the maximum modulus of particle
displacements.

*/
class NearestNeighbours {
public:
  struct Pair {int first,second; Pair(int i,int j) : first(i),second(j) {} } ;

  NearestNeighbours(double rc_) : topological(false), rc(rc_) {}
  NearestNeighbours(int Nnearest_) : topological(true), Nnearest(Nnearest_) {}
  virtual ~NearestNeighbours() {}
  virtual void rebuild(OLconfiguration&,double rc=-1)=0;  ///< Build lists from scratch for given conf
  virtual void rebuild(OLconfiguration&,int n)=0; ///< Build list from scratch for n neighbour case, does not build pair list since it has no meaning in this case
  virtual void update(OLconfiguration&,double)=0; ///< Inform of change in configuration, will try to update lists assuming particles have not moved much, may rebuild everything from scratch

  virtual std::vector<Pair>::iterator pair_begin()=0;
  virtual std::vector<Pair>::iterator pair_end()=0;
  virtual std::vector<int>::iterator neighbours_begin(int i)=0;
  virtual std::vector<int>::iterator neighbours_end(int i)=0;

protected:
  bool   topological;
  int    Nnearest,
  double rc;
} ;


/** \class NeighbourList_naive
    \ingroup OfflatticeINT

    This provides a naive implementation of a neighbour list.
    Although the list is built by iterating through all possible
    pairs, if the simulation step provides information on the maximum
    displacement, this brings noticeable speed improvements at more
    than 200 particles, at least for potentials with relatively short
    cut-offs (such as the repulsive Lennard-Jones).
*/
class NeighbourList_naive : public NearestNeighbours {
public:
  NeighbourList_naive(double rc,double delta_r=-1);
  void rebuild(OLconfiguration&,double rc=-1);
  void update(OLconfiguration&,double maxdisp);

  std::vector<Pair>::iterator pair_begin() {return pairs.begin();}
  std::vector<Pair>::iterator pair_end() {return pairs.end();}
  std::vector<int>::iterator neighbours_begin(int i) {return neighbours[i].begin();}
  std::vector<int>::iterator neighbours_end(int i) {return neighbours[i].end();}

private:
  std::vector<Pair> pairs;
  std::vector<std::vector<int>>  neighbours;

  double rdsq,delta_r;
  double accum_maxdisp;
} ;

/** \class TopologicalNeighbours_naive
    \ingroup OfflatticeINT

    This produces pairs of nearests neighbours, found with the usual
    Euclidean distance, but for each particle finds exactly N nearest
    neighbours, thus leading to potentially non-symmetric interaction
    matrices.  Implementation is naive (linear search of neighbours
    for each particle).  It is rather slow and meant mainly as a check
    for implementations with better algorithms.
*/
class TopologicalNeighbours_naive : public NearestNeighbours {
public:
  TopologicalNeighbours_naive(int Nnearest) :
    NearestNeighbours(Nnearest) {};
  void rebuild(OLconfiguration&,int Nnearest=-1);
  void update(OLconfiguration&,double maxdisp);

  std::vector<Pair>::iterator pair_begin() {return pairs.begin();}
  std::vector<Pair>::iterator pair_end() {return pairs.end();}
  std::vector<int>::iterator neighbours_begin(int i) {return neighbours[i].begin();}
  std::vector<int>::iterator neighbours_end(int i) {return neighbours[i].end();}

private:
  struct ndist {int n; double dsq;};
  std::vector<Pair> pairs;
  //  std::vector<std::vector<int>>  neighbours;
  std::vector<std::vector<int>>  neighbours;

  int Nnearest;
} ;


} /* namespace */

#endif /* NNEIGHBOURS_HH */
