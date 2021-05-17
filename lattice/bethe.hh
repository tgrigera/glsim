/*
 * bethe.hh -- classes for the Bethe lattice
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

#ifndef BETHE_HH
#define BETHE_HH

#include "graph.hh"

namespace glsim {

/** \ingroup lattice
 * \brief class for the Bethe lattice
 *
 * The Bethe lattice is a graph with fixed connectivity but no loops
 * (i.e. it is a tree where the root has z children, and each nonroot
 * node has z-1 children).  In a finite realization, as here, each
 * level of the tree has z neighbours except the lowest, which has
 * only one.
 *
 * A Bethe lattice with \f$ L\f$ layers (appart from the central node,
 * or layer 0) has 1 node at layer 0, \f$z\f$ nodes at layer 1,
 * \f$z(z-1)\f$ nodes at layer 2.  In general, for \f$k\ge1\f$, layer
 * \f$k\f$ has \f$ N_k = z(z-1)^k \f$ nodes.  The total number of
 * nodes is thus \f[ N = 1 + z \frac{(z-1)^L-1}{z-2}. \f]
 *
 */
template <typename nodeT>
class BetheLattice : public GraphBase<nodeT> {
public:
  BetheLattice(int z,int L);  // L is number of layers (not counting
			      // the central node which is layer 0)
private:
  int  coordination;
} ;

/**
 * \param z is the coordination number
 * \param L is the number of layers (not counting the base or central node)
 */
template <typename nodeT>
BetheLattice<nodeT>::BetheLattice(int z,int L) :
  GraphBase<nodeT>(),
  coordination(z)
{
  int N=1;
  int Nlast_layer=z;
  for (int k=1; k<=L; ++k) {
    N+=Nlast_layer;
    Nlast_layer*=z-1;
  }
  Nlast_layer/=z-1;
  GraphBase<nodeT>::init_storage(N,-1,true);

  id_t last;
  for (last=1; last<=z; ++last)
    GraphBase<nodeT>::add_bond(0,last);
  for (id_t i=1; i<N-Nlast_layer; ++i)
    for (int nn=1; nn<=z-1; ++nn)
      GraphBase<nodeT>::add_bond(i,last++);
}

} /* namespace */

#endif /* BETHE_HH */
