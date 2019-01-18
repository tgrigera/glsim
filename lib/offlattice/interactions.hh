/*
 * interactions.hh -- Interactions for off-lattice systems
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

#ifndef INTERACTIONS_HH
#define INTERACTIONS_HH

#include "olconfiguration.hh"

#include "nneighbours.hh"

namespace glsim {
  
/*****************************************************************************
 *
 * Interactions
 *
 */

/** \class Interactions
    \ingroup OfflatticeINT

The interactions object here is abstract.  We allow the constructor to
take a configuration as argument, so that consistency checks on the
configuration can be done.  Thus Interactions must be created _after_
loading configuration.  But there will be a Potential object to be
created before reading environment (since it can create its own
private environment to be initialized by prepare).


To write a potential object

 - Provide functions xxx

 - Create potential before creating interactions.

 - Potential need not be completely initialized; interactions will
   call a potential init that may also do checks on configuration.
   This way the potential constructor can create its own environment.

*/
class Interactions {
public:
  Interactions() {}
  virtual const char *name() const =0;
  virtual bool   conserve_P() const =0; ///< Whether the interactions conserves momentum (e.g. has external field)
  virtual double mass(short type) const =0;
  virtual double potential_energy(OLconfiguration&) = 0;
  virtual double force_and_potential_energy(OLconfiguration&) = 0;
  virtual double acceleration_and_potential_energy(OLconfiguration&) = 0;
  virtual double kinetic_energy_and_momentum(OLconfiguration&,double P[]);
  virtual double delta_energy_particle_shift(OLconfiguration&,int,double*) = 0;
  virtual double delta_energy_particle_swap(OLconfiguration&,int a,int b) = 0;

  virtual void   fold_coordinates(OLconfiguration&,double maxdisp=-1.);
  virtual ~Interactions() {}
} ;
  
/// After moving particles, call this.  maxdisp is a measure of how
/// much particles have moved in this step (mostly an upper bound on
/// displacemets).  Interactions keeps track of these movements and
/// reinits structures as needed.  If -1 means reset now.
inline void Interactions::fold_coordinates(OLconfiguration& conf,double maxdisp)
{
  conf.fold_coordinates();
}

/** \class FreeParticles
    \ingroup OfflatticeINT

    Trivial implementation of interactions for free particles.
*/
class FreeParticles : public Interactions {
  const char* name() const {return "Free particles";}
  double mass(short t) const {return 1.;}
  bool   conserve_P() const {return true;}
  double potential_energy(OLconfiguration& c) {return 0;}
  double force_and_potential_energy(OLconfiguration& c) {return 0;}
  double acceleration_and_potential_energy(OLconfiguration& c) {return 0;}
  double delta_energy_particle_shift(OLconfiguration&,int,double*) {return 0;}
  double delta_energy_particle_swap(OLconfiguration&,int,int) {return 0;}
} ;

/*****************************************************************************
 *
 * Interactions_pairwise_isotropic
 *
 */

/** \class Interactions_isotropic_pairwise_naive.
    \ingroup OfflatticeINT

This class can compute force and energy for isotropic pairwise
interactions.  We will provide several versions of this, since there
are several methods to find all pairs within the cut-off radius
without going through all \f$ N(N-1)/2 \f$ possibilities. We start with a
version that always goes through all the pairs. It is the easiest to
write and will be useful as a minimal implementation against which
more intelligent algorithms can be tested, and it will be actually
used in simulations of very small systems.

    We write these as templates so that inlined functions can be used when
inserting the actual pair potential.

*/
template <typename PotentialT>
class Interactions_isotropic_pairwise_naive : public Interactions {
public:
  Interactions_isotropic_pairwise_naive(PotentialT &P,OLconfiguration&);
  const  char *name() const;
  double mass(short type) const {return PP.mass(type);}
  bool   conserve_P() const {return !PP.has_efield();}
  double potential_energy(OLconfiguration&);
  double force_and_potential_energy(OLconfiguration&);
  double acceleration_and_potential_energy(OLconfiguration&);
  double delta_energy_particle_shift(OLconfiguration &conf,int n,double *rnew);
  double delta_energy_particle_swap(OLconfiguration&,int a,int b);
  void   tabulate_potential(std::ostream&,short t1,short t2);

private:
  PotentialT& PP;
} ;

template <typename PotentialT> inline
Interactions_isotropic_pairwise_naive<PotentialT>::
Interactions_isotropic_pairwise_naive(PotentialT &p,OLconfiguration& c) :
  Interactions(), PP(p)
{
  p.init(c);
}

template <typename PotentialT> inline
const char *Interactions_isotropic_pairwise_naive<PotentialT>::name() const
{
  return PP.name();
}

template <typename PotentialT>
double Interactions_isotropic_pairwise_naive<PotentialT>::
potential_energy(OLconfiguration& conf)
{
  int    n,m;
  double E;

  E=0;
  for (n=0; n<conf.N; n++) {
    double rn[3];
    rn[0]=conf.r[n][0];
    rn[1]=conf.r[n][1];
    rn[2]=conf.r[n][2];
    int    typen=conf.type[n];

    if (PP.has_efield())
      E+=PP.external_field(rn,typen);
    
    for (m=n+1; m<conf.N; m++) {
      int    typem=conf.type[m];
      double dsq=conf.distancesq(rn,conf.r[m]);
      if (!PP.within_cutoff(dsq,typen,typem)) continue;
      E+=PP.pair_potential_r(dsq,typen,typem);
    }
  }
  
  return E;
}

/**
   NOTE that we expect v'(r)/r rather than the derivative

*/
template <typename PotentialT>
double Interactions_isotropic_pairwise_naive<PotentialT>::
force_and_potential_energy(OLconfiguration &conf)
{
  int    n,m;
  double vpr,Et,E=0,force[3];
  memset(conf.a,0,3*conf.N*sizeof(double));

  for (int n=0; n<conf.N; n++) { 
    double rn[3];
    memcpy(rn,conf.r[n],3*sizeof(double));
    int typen=conf.type[n];

    if (PP.has_efield()) {
      E+=PP.external_field(rn,typen,force);
      conf.a[n][0]+=force[0];
      conf.a[n][1]+=force[1];
      conf.a[n][2]+=force[2];
    }

    for (int m=n+1; m<conf.N; m++) {
      int    typem=conf.type[m];
      double rxmn=conf.ddiff(conf.r[m][0],rn[0],conf.box_length[0]);
      double rymn=conf.ddiff(conf.r[m][1],rn[1],conf.box_length[1]);
      double rzmn=conf.ddiff(conf.r[m][2],rn[2],conf.box_length[2]);
      double dsq=rxmn*rxmn+rymn*rymn+rzmn*rzmn;
      if (!PP.within_cutoff(dsq,typen,typem)) continue;
      E+=PP.pair_potential_r(dsq,typen,typem,vpr);
      conf.a[m][0]-=vpr*rxmn;
      conf.a[m][1]-=vpr*rymn;
      conf.a[m][2]-=vpr*rzmn;
      conf.a[n][0]+=vpr*rxmn;
      conf.a[n][1]+=vpr*rymn;
      conf.a[n][2]+=vpr*rzmn;
    }
  }

  return E;
}

template <typename PotentialT>
double Interactions_isotropic_pairwise_naive<PotentialT>::
acceleration_and_potential_energy(OLconfiguration &conf)
{
  int    n,m;
  double force[3],vpr,Et,E=0;
  memset(conf.a,0,3*conf.N*sizeof(double));

  for (int n=0; n<conf.N-1; n++) { 
    double rn[3];
    memcpy(rn,conf.r[n],3*sizeof(double));
    int typen=conf.type[n];

    if ( (Et=PP.external_field(rn,typen,force))!=0) {
      E+=Et;
      conf.a[n][0]+=force[0]/PP.mass(typen);
      conf.a[n][1]+=force[1]/PP.mass(typen);
      conf.a[n][2]+=force[2]/PP.mass(typen);
    }

    for (int m=n+1; m<conf.N; m++) {
      int    typem=conf.type[m];
      double rxmn=conf.ddiff(conf.r[m][0],rn[0],conf.box_length[0]);
      double rymn=conf.ddiff(conf.r[m][1],rn[1],conf.box_length[1]);
      double rzmn=conf.ddiff(conf.r[m][2],rn[2],conf.box_length[2]);
      double dsq=rxmn*rxmn+rymn*rymn+rzmn*rzmn;
      if (!PP.within_cutoff(dsq,typen,typem)) continue;
      E+=PP.pair_potential_r(dsq,typen,typem,vpr);
      double vprm=vpr/PP.mass(typem);
      double vprn=vpr/PP.mass(typen);
      conf.a[m][0]-=vprm*rxmn;
      conf.a[m][1]-=vprm*rymn;
      conf.a[m][2]-=vprm*rzmn;
      conf.a[n][0]+=vprn*rxmn;
      conf.a[n][1]+=vprn*rymn;
      conf.a[n][2]+=vprn*rzmn;
    }
  }

  return E;
}

template <typename PotentialT>
double Interactions_isotropic_pairwise_naive<PotentialT>::
delta_energy_particle_shift(OLconfiguration &conf,int n,double *rnew)
{
  double rn[3];
  int typen=conf.type[n];
  memcpy(rn,conf.r[n],3*sizeof(double));

  double PE=PP.external_field(rn,typen);
  double PEshift=PP.external_field(rnew,typen);
  double DE=PEshift-PE;
  for (int m=0; m<conf.N; m++) {
    if (m==n) continue;
    int    typem=conf.type[m];
    double dsq=conf.distancesq(rn,conf.r[m]);
    PE=PP.pair_potential(dsq,typen,typem);
    dsq=conf.distancesq(rnew,conf.r[m]);
    PEshift=PP.pair_potential(dsq,typen,typem);
    DE+=PEshift-PE;
  }
  return DE;
}

template <typename PotentialT>
double Interactions_isotropic_pairwise_naive<PotentialT>::
delta_energy_particle_swap(OLconfiguration& conf,int a,int b)
{
  double rn[3],E0=0,E1=0;

  int typea=conf.type[a];
  int typeb=conf.type[b];

  if (PP.has_efield()) {
    E0=PP.external_field(conf.r[a],typea);
    E0+=PP.external_field(conf.r[b],typeb);
    E1=PP.external_field(conf.r[a],typeb);
    E1+=PP.external_field(conf.r[b],typea);
  }

  for (int m=0; m<conf.N; m++) {
    if (m==a || m==b) continue;    // no int with itself; no difference in DE
                                  // by omitting the pair a-b from the loop
    int    typem=conf.type[m];
    double dsq=conf.distancesq(a,m);
    E0+=PP.pair_potential(dsq,typea,typem);
    E1+=PP.pair_potential(dsq,typeb,typem);
  }
  for (int m=0; m<conf.N; m++) {
    if (m==a || m==b) continue;
    int    typem=conf.type[m];
    double dsq=conf.distancesq(b,m);
    E0+=PP.pair_potential(dsq,typeb,typem);
    E1+=PP.pair_potential(dsq,typea,typem);
  }

  return E1-E0;
}  




template <typename PotentialT>
void Interactions_isotropic_pairwise_naive<PotentialT>::tabulate_potential(std::ostream& o,
									  short t1,short t2)
{
  int    N=200;
  double r1=PP.cutoff();
  double delta=r1/N;
  double r0=delta;

  o << "# r    V(r) for types " << t1 << " and " << t2 << "\n#\n";
  for (int i=0; i<N+10; i++) {
    double r=r0+delta*i;
    o << r << "  " << PP.pair_potential(r*r,t1,t2) << '\n';
  }
}

/*****************************************************************************
 *
 * class Interactions_isotropic_pairwise
 *
 */

/** \class  Interactions_isotropic_pairwise
    \ingroup OfflatticeINT

    This class provides functions to compute energy, acceleration and
force given an isotropic pair potential that must be indicated as a
template parameter.  It uses an object of type NearestNeighbours to
find the particles within the cut-off distance.  

*/
template <typename PotentialT,typename NeighboursT=NeighbourList_subcells>
class Interactions_isotropic_pairwise : public Interactions {
public:
  Interactions_isotropic_pairwise(PotentialT &P,OLconfiguration &c,
				  NeighboursT *NN=0);
  ~Interactions_isotropic_pairwise();
  const  char *name() const;
  double mass(short type) const {return PP.mass(type);}
  bool   conserve_P() const {return !PP.has_efield();}
  double potential_energy(OLconfiguration&);
  double force_and_potential_energy(OLconfiguration&);
  double acceleration_and_potential_energy(OLconfiguration&);
  double delta_energy_particle_shift(OLconfiguration &conf,int n,double *rnew);
  double delta_energy_particle_swap(OLconfiguration &conf,int n,int m);
  void   tabulate_potential(std::ostream&,short t1,short t2);
  void   fold_coordinates(OLconfiguration& conf,double maxdisp=-1);

private:
  PotentialT&  PP;
  bool         own_NN;
  NeighboursT* NN;
} ;

template <typename PotentialT,typename NeighboursT> inline
void Interactions_isotropic_pairwise<PotentialT,NeighboursT>::
fold_coordinates(OLconfiguration& conf,double maxdisp)
{
  conf.fold_coordinates();
  if (maxdisp<0)
    NN->rebuild(conf,PP.cutoff());
  else
    NN->update(maxdisp);
}

/**
The last argument is a pointer to an Neighbour object (see
nneighbours.hh).  If null (default), an object will be created passing
only the cutoff as argument.
*/
template <typename PotentialT,typename NeighboursT> inline
Interactions_isotropic_pairwise<PotentialT,NeighboursT>::
Interactions_isotropic_pairwise(PotentialT &p,OLconfiguration& c,
				NeighboursT *n) :
  Interactions(), PP(p)
{
  PP.init(c);
  if (n) {
    NN=n;
    own_NN=false;
  } else {
    NN=new NeighboursT(PP.cutoff());
    own_NN=true;
  }
  NN->rebuild(c,PP.cutoff());
}

template <typename PotentialT,typename NeighboursT> inline
Interactions_isotropic_pairwise<PotentialT,NeighboursT>::
~Interactions_isotropic_pairwise()
{
  if (own_NN) delete NN;
}

template <typename PotentialT,typename NeighboursT> inline
const char *Interactions_isotropic_pairwise<PotentialT,NeighboursT>::name() const
{
  return PP.name();
}

template <typename PotentialT,typename NeighboursT>
double Interactions_isotropic_pairwise<PotentialT,NeighboursT>::
potential_energy(OLconfiguration& conf)
{
  int    n,m;
  double E;

  E=0;
  if (PP.has_efield())
    for (n=0; n<conf.N; n++)
      E+=PP.external_field(conf.r[n],conf.type[n]);

  for (auto p = NN->pairs_begin(); p!=NN->pairs_end(); ++p) {
      int    typen=conf.type[p->first];
      int    typem=conf.type[p->second];
      double dsq=conf.distancesq(conf.r[p->first],conf.r[p->second]);
      if (!PP.within_cutoff(dsq,typen,typem)) continue;
      E+=PP.pair_potential_r(dsq,typen,typem);
  }
  
  return E;
}

/**
   NOTE that we expect v'(r)/r rather than the derivative

*/
template <typename PotentialT,typename NeighboursT>
double Interactions_isotropic_pairwise<PotentialT,NeighboursT>::
force_and_potential_energy(OLconfiguration &conf)
{
  int    n,m;
  double vpr,Et,E=0,force[3];
  memset(conf.a,0,3*conf.N*sizeof(double));

  if (PP.has_efield()) {
    for (int n=0; n<conf.N; n++) { 
      E+=PP.external_field(conf.r[n],conf.type[n],force);
      conf.a[n][0]+=force[0];
      conf.a[n][1]+=force[1];
      conf.a[n][2]+=force[2];
    }
  }

  for (auto p = NN->pairs_begin(); p!=NN->pairs_end(); ++p) {
    int typen=conf.type[p->first];
    int typem=conf.type[p->second];
    double rxmn=conf.ddiff(conf.r[p->first][0],conf.r[p->second][0],conf.box_length[0]);
    double rymn=conf.ddiff(conf.r[p->first][1],conf.r[p->second][1],conf.box_length[1]);
    double rzmn=conf.ddiff(conf.r[p->first][2],conf.r[p->second][2],conf.box_length[2]);
    double dsq=rxmn*rxmn+rymn*rymn+rzmn*rzmn;
    if (!PP.within_cutoff(dsq,typen,typem)) continue;
    E+=PP.pair_potential_r(dsq,typen,typem,vpr);
    conf.a[p->first][0]-=vpr*rxmn;
    conf.a[p->first][1]-=vpr*rymn;
    conf.a[p->first][2]-=vpr*rzmn;
    conf.a[p->second][0]+=vpr*rxmn;
    conf.a[p->second][1]+=vpr*rymn;
    conf.a[p->second][2]+=vpr*rzmn;
  }

  return E;
}

template <typename PotentialT,typename NeighboursT>
double Interactions_isotropic_pairwise<PotentialT,NeighboursT>::
acceleration_and_potential_energy(OLconfiguration &conf)
{
  int    n,m;
  double vpr,Et,E=0,force[3];
  memset(conf.a,0,3*conf.N*sizeof(double));

  if (PP.has_efield()) {
    for (int n=0; n<conf.N; n++) { 
      E+=PP.external_field(conf.r[n],conf.type[n],force);
      conf.a[n][0]+=force[0]/PP.mass(conf.type[n]);
      conf.a[n][1]+=force[1]/PP.mass(conf.type[n]);
      conf.a[n][2]+=force[2]/PP.mass(conf.type[n]);
    }
  }

  for (auto p = NN->pairs_begin(); p!=NN->pairs_end(); ++p) {
    int typen=conf.type[p->first];
    int typem=conf.type[p->second];
    double rxmn=conf.ddiff(conf.r[p->first][0],conf.r[p->second][0],conf.box_length[0]);
    double rymn=conf.ddiff(conf.r[p->first][1],conf.r[p->second][1],conf.box_length[1]);
    double rzmn=conf.ddiff(conf.r[p->first][2],conf.r[p->second][2],conf.box_length[2]);
    double dsq=rxmn*rxmn+rymn*rymn+rzmn*rzmn;
    if (!PP.within_cutoff(dsq,typen,typem)) continue;
    E+=PP.pair_potential_r(dsq,typen,typem,vpr);
    double vprm=vpr/PP.mass(typem);
    double vprn=vpr/PP.mass(typen);
    conf.a[p->first][0]-=vpr*rxmn;
    conf.a[p->first][1]-=vpr*rymn;
    conf.a[p->first][2]-=vpr*rzmn;
    conf.a[p->second][0]+=vpr*rxmn;
    conf.a[p->second][1]+=vpr*rymn;
    conf.a[p->second][2]+=vpr*rzmn;
  }

  return E;
}

template <typename PotentialT,typename NeighboursT>
double Interactions_isotropic_pairwise<PotentialT,NeighboursT>::
delta_energy_particle_shift(OLconfiguration &conf,int n,double *rnew)
{
  double PE,PEshift,rn[3];
  int typen=conf.type[n];
  memcpy(rn,conf.r[n],3*sizeof(double));

  double DE=0;
  if (PP.has_efield()) {
    PE=PP.external_field(rn,typen);
    PEshift=PP.external_field(rnew,typen);
    DE=PEshift-PE;
  }
  for (auto m=NN->neighbours_begin(n); m!=NN->neighbours_end(n); ++m) {
    int    typem=conf.type[*m];
    double dsq=conf.distancesq(rn,conf.r[*m]);
    PE=PP.pair_potential(dsq,typen,typem);
    dsq=conf.distancesq(rnew,conf.r[*m]);
    PEshift=PP.pair_potential(dsq,typen,typem);
    DE+=PEshift-PE;
  }
  return DE;
}


template <typename PotentialT,typename NeighboursT>
double Interactions_isotropic_pairwise<PotentialT,NeighboursT>::
delta_energy_particle_swap(OLconfiguration &conf,int a,int b)
{
  double E0=0,E1=0;
  int typea=conf.type[a];
  int typeb=conf.type[b];

  if (PP.has_efield()) {
    E0=PP.external_field(conf.r[a],typea);
    E0+=PP.external_field(conf.r[b],typeb);
    E1=PP.external_field(conf.r[a],typeb);
    E1+=PP.external_field(conf.r[b],typea);
  }

  for (auto m=NN->neighbours_begin(a); m!=NN->neighbours_end(a); ++m) {
    int    typem=conf.type[*m];
    double dsq=conf.distancesq(a,*m);
    E0+=PP.pair_potential(dsq,typea,typem);
    E1+=PP.pair_potential(dsq,typeb,typem);
  }
  for (auto m=NN->neighbours_begin(b); m!=NN->neighbours_end(b); ++m) {
    int    typem=conf.type[*m];
    double dsq=conf.distancesq(b,*m);
    E0+=PP.pair_potential(dsq,typeb,typem);
    E1+=PP.pair_potential(dsq,typea,typem);
  }

  return E1-E0;
  
}


} /* namespace */

#endif /* INTERACTIONS_HH */
