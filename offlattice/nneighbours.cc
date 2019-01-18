/*
 * nnneighbours.cc -- structures to find nearest neighbours
 *
 * This file is part of olglsim, a numerical simulation class library
 * and helper programs.
 *
 * olglsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * olglsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use olglsim to produced published work, or if you redistribute a
 * modified version of olglsim, or code based on olglsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * olglsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding.
 *
 * olglsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

#include <queue>
#include "glsim/log.hh"
#include "nneighbours.hh"

namespace glsim {

/*****************************************************************************
 *
 * NeighbourList_naive
 *
 */

NeighbourList_naive::NeighbourList_naive(double rc_,double delta_r_) :
  rc(rc_),
  delta_r(delta_r_),
  accum_maxdisp(0.)
{
  if (delta_r<0)
    delta_r=rc*0.3;
  rcsq=rc*rc;
  rdsq=rc+delta_r;
  rdsq*=rdsq;
}

void NeighbourList_naive::rebuild(OLconfiguration& conf_,double rc_)
{
  conf=&conf_;
  if (rc_>0) {
    rc=rc_;
    rcsq=rc*rc;
    rdsq=rc+delta_r;
    rdsq*=rdsq;
  }
  
  neighbours.clear();
  neighbours.resize(conf->N);
  pairs.clear();
  for (size_t i=0; i<conf->N-1; ++i)
    for (size_t j=i+1; j<conf->N; ++j)
      if (conf->distancesq(i,j) <= rdsq) {
	pairs.push_back(Pair(i,j));
	neighbours[i].push_back(j);
	neighbours[j].push_back(i);
      }
  accum_maxdisp=0.;
}

void NeighbourList_naive::update(double maxdisp)
{
  accum_maxdisp+=maxdisp;
  if (accum_maxdisp>=delta_r/2.)
    rebuild(*conf);
}


/******************************************************************************
 *
 * Subcells (metric nearest neighbours with subcell method)
 *
 */

void Subcells::init_cells(int npart,double box_length[])
{
  nparticles=npart;
  /* Initialization of subcell index */
  /* The factor 1.001 is a "security factor", to provide for the case when
     numerical precision makes a particle stay just outside the periodic box */
  float sfac=1.001;
  double scmin=rc+delta_r;
  for (int i=0; i<3; i++) {
    m[i]=(int) floor(box_length[i]/scmin);
    boxl[i]=box_length[i];
    if (m[i]<3) {
      std::cerr << "Got " << m[i] << " subcells!\n";
      throw System_too_small(i==0 ? "x" : (i==1 ? "y" : "z"),boxl[i],3*scmin);
    } else scell[i]=sfac*box_length[i]/m[i]; 
  }
  nscell=m[0]*m[1]*m[2];
  glsim::logs(glsim::info) << "Subcells for nearest neighbours reset, using " <<
    m[0] << " x " << m[1] << " x " << m[2] << " cells.\n";

  clear_lists();
  subcell=new int[nscell+1];  // last subcell is an end marker (one-past-end)
  subcell[nscell]=-1;
  for (int i=0; i<nscell; i++) subcell[i]=-1;

  llist=new int[nparticles];
  WhichSubcell=new int[nparticles];
  subcelln=new int[27*nscell];

  /* Initialization of subcell neighbour list */
  for (int ix=0; ix<m[0]; ix++)
    for (int iy=0; iy<m[1]; iy++)
      for (int iz=0; iz<m[2]; iz++) {
	int i=icell(ix,iy,iz)*27;
	subcelln[i]   = eicell(ix-1,iy+1,iz  );
	subcelln[i+1] = eicell(ix  ,iy+1,iz  );
	subcelln[i+2] = eicell(ix+1,iy+1,iz  );
	subcelln[i+3] = eicell(ix+1,iy  ,iz  );
	subcelln[i+4] = eicell(ix-1,iy+1,iz-1);
	subcelln[i+5] = eicell(ix  ,iy+1,iz-1);
	subcelln[i+6] = eicell(ix+1,iy+1,iz-1);
	subcelln[i+7] = eicell(ix-1,iy  ,iz-1);
	subcelln[i+8] = eicell(ix  ,iy  ,iz-1);
	subcelln[i+9] = eicell(ix+1,iy  ,iz-1);
	subcelln[i+10]= eicell(ix-1,iy-1,iz-1);
	subcelln[i+11]= eicell(ix  ,iy-1,iz-1);
	subcelln[i+12]= eicell(ix+1,iy-1,iz-1);
        /* Neighbours 13 and up are only used in Monte Carlo or for
           border subcells */
	subcelln[i+13]= eicell(ix-1,iy-1,iz  );
	subcelln[i+14]= eicell(ix  ,iy-1,iz  );
	subcelln[i+15]= eicell(ix+1,iy-1,iz  );
	subcelln[i+16]= eicell(ix-1,iy  ,iz  );
	subcelln[i+17]= eicell(ix-1,iy+1,iz+1);
	subcelln[i+18]= eicell(ix  ,iy+1,iz+1);
	subcelln[i+19]= eicell(ix+1,iy+1,iz+1);
	subcelln[i+20]= eicell(ix-1,iy  ,iz+1);
	subcelln[i+21]= eicell(ix  ,iy  ,iz+1);
	subcelln[i+22]= eicell(ix+1,iy  ,iz+1);
	subcelln[i+23]= eicell(ix-1,iy-1,iz+1);
	subcelln[i+24]= eicell(ix  ,iy-1,iz+1);
	subcelln[i+25]= eicell(ix+1,iy-1,iz+1);
	subcelln[i+26]= nscell;  // This is a one-past-end subcell for convinence in operator++
  }
}

void Subcells::clear_lists()
{
  delete[] subcell; subcell=0;
  delete[] llist; llist=0;
  delete[] WhichSubcell; WhichSubcell=0;
  delete[] subcelln; subcelln=0;
}

void Subcells::rebuild(OLconfiguration& conf_,double rc_,double delta_rc_)
{
  conf=&conf_;
  bool rcchange = rc_>0 && rc_!=rc;
  if (rcchange) {
    rc=rc_;
    rcsq=rc*rc;
  }
  rcchange = rcchange || (delta_rc_>0 && delta_rc_!=delta_r);
  if (delta_rc_>0) delta_r=delta_rc_;

  if (rcchange || nparticles!=conf->N || conf->box_length[0]!=boxl[0] ||
      conf->box_length[1]!=boxl[1] || conf->box_length[2]!=boxl[2])
    init_cells(conf->N,conf->box_length);

  update(delta_r);
}

void Subcells::update(double maxdisp)
{
  int i,ix,iy,iz,ic;

  accum_maxdisp+=maxdisp;
  if (accum_maxdisp<delta_r/2.) return;
  accum_maxdisp=0;

  for (i=0; i<nscell; i++)
    subcell[i]=-1;
  for (i=0; i<conf->N; i++) {
    ix=(int) floor(conf->r[i][0]/scell[0]);
    iy=(int) floor(conf->r[i][1]/scell[1]);
    iz=(int) floor(conf->r[i][2]/scell[2]);
    ic=icell(ix,iy,iz);

#ifdef DEBUG
    if (ix<0 || ix>=m[0] || iy<0 || iy>=m[1] || iz<0 || iz>=m[2]) {
      std::cout << "ERROR in scell: ix iy iz " << ix << ' ' << iy << ' ' << iz << '\n';
      std::cout << "                mx my mz " << m[0] << ' ' << m[1] << ' ' << m[2] << '\n';
      std::cout << "                particle " << i << '\n';
      std::cout << "                position " << conf->r[i][0] << ' ' << conf->r[i][1] << ' ' << conf->r[i][2] << '\n';
      throw glsim::Runtime_error("Invalid subcell number",HERE);
    }
#endif

    llist[i]=subcell[ic];
    subcell[ic]=i;
    WhichSubcell[i]=ic;
  }
}

Subcells::NeighbourIterator& Subcells::NeighbourIterator::operator++()
{
  n=scell.next_particle(n);
  if (n==particle) n=scell.next_particle(n);
  while (n<0 && ncell<26) {
    ncell++;
    n=scell.first_particle(scell.neighbour(icell,ncell));
  }
  return *this;
};

Subcells::PairIterator& Subcells::PairIterator::operator++()
{
  pair.second=scell.next_particle(pair.second);
  while (pair.second<0 && ncell<13) {
    ncell++;
    pair.second=scell.first_particle(scell.neighbour(icell,ncell));
  }
  if (pair.second<0 || ncell>=13) {  // Neighbours of iparticle exhausted, get next iparticle
    ncell=-1;
    pair.first=scell.next_particle(pair.first);
    while (pair.first<0 && icell<scell.num_subcells() ) {
      icell++;
      pair.first=scell.first_particle(icell);
    }
    if (pair.first>=0) {  // If iparticles not exhausted, get next jparticle
      pair.second=scell.next_particle(pair.first);
      while (pair.second<0 && ncell<13) {
	ncell++;
	pair.second=scell.first_particle(scell.neighbour(icell,ncell));
      }
      if (pair.second<0) pair.first=-1;  // Reached the end
    } else pair.second=-1; // Reached the end
  }
  return *this;
};


/*****************************************************************************
 *
 * NeighbourList_subcells
 *
 */

NeighbourList_subcells::NeighbourList_subcells(double rc_,double delta_r_) :
  rc(rc_),
  delta_r(delta_r_),
  accum_maxdisp(0.),
  SUBC(rc_+delta_r_,0)
{
  if (delta_r<0)
    delta_r=rc*0.3;
  rcsq=rc*rc;
  rdsq=rc+delta_r;
  rdsq*=rdsq;
}

void NeighbourList_subcells::rebuild(OLconfiguration& conf_,double rc_)
{
  conf=&conf_;
  if (rc_>0) {
    rc=rc_;
    rcsq=rc*rc;
    rdsq=rc+delta_r;
    rdsq*=rdsq;
  }
  
  neighbours.clear();
  neighbours.resize(conf->N);
  pairs.clear();
  SUBC.rebuild(*conf,rc+delta_r,0);
  for_each_pair(SUBC,
		[this](int i,int j,double dsq){
		  pairs.push_back(Pair(i,j));
		  neighbours[i].push_back(j);
		  neighbours[j].push_back(i);
		}
		);
  accum_maxdisp=0.;
}

void NeighbourList_subcells::update(double maxdisp)
{
  accum_maxdisp+=maxdisp;
  if (accum_maxdisp>=delta_r/2.)
    rebuild(*conf);
}


/*****************************************************************************
 *
 * TopologicalNeighbours_naive
 *
 */

void TopologicalNeighbours_naive::rebuild(OLconfiguration& conf,int Nnearest_)
{
  if (Nnearest_>0) Nnearest=Nnearest_;

  struct comp {
    bool operator()(ndist& a,ndist& b) {return a.dsq<b.dsq;}
  } ;

  neighbours.clear();
  neighbours.resize(conf.N);

  std::priority_queue<ndist,std::vector<ndist>,comp>  candidates;

  for (size_t i=0; i<conf.N; ++i) {
    for (size_t j=0; j<conf.N; ++j) {
      if (i==j) continue;
      double dsq=conf.distancesq(i,j);
      if (candidates.size()<Nnearest)
	candidates.push(ndist(j,dsq));
      else if (dsq<candidates.top().dsq) {
	candidates.pop();
	candidates.push(ndist(j,dsq));
      }
    }
    while (!candidates.empty()) {
      neighbours[i].push_back(candidates.top().n);
      candidates.pop();
    }
  }
}


#ifdef UFA
class TopologicalNeighbours_naive::neighbour_iterator :
  public std::iterator<std::forward_iterator_tag,int>
{
public:
  ///\name Construction, copy and comparison
  //@{
  nieghbour_iterator(PeriodicSQLattice<nodeT> &l,id_t i_=0,id_t j_=0);
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

#endif













} /* namespace */
