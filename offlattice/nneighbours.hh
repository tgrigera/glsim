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

class NearestNeighbours {
public:
  struct Pair {int first,second; Pair(int i,int j) : first(i),second(j) {} } ;

  NearestNeighbours() {}
  void reset(OLconfiguration&,double rsq);

  std::vector<Pair>::iterator pair_begin() {return pairs.begin();}
  std::vector<Pair>::iterator pair_end() {return pairs.end();}
  std::vector<int>::iterator neighbours_begin(int i) {return neighbours[i].begin();}
  std::vector<int>::iterator neighbours_end(int i) {return neighbours[i].end();}

private:
  std::vector<Pair> pairs;
  std::vector<std::vector<int>>  neighbours;
} ;


} /* namespace */

#endif /* NNEIGHBOURS_HH */
