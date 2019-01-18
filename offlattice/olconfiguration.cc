/*
 * olconfiguration.cc -- definitions for class OLconfiguration
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

#include <string.h>
#include <cmath>
#include <utility>
#include <unordered_map>

#include "glsim/cerrors.h"
#include "olconfiguration.hh"

namespace glsim {

/*****************************************************************************
 *
 * OLconfiguration
 *
 */


/*
 * Construction, destruction, copy and swap
 *
 */

OLconfiguration::~OLconfiguration()
{
  delete[] id;
  delete[] type;
  delete[] flags;
  delete[] r;
  delete[] v;
  delete[] a;
  delete[] reference_r;
  N=0;
}

OLconfiguration::OLconfiguration(const OLconfiguration& o) :
  Configuration(o.name),
  N(o.N),
  time(o.time),
  step(o.step),
  id(0), type(0), flags(0), r(0), v(0), a(0),
  reference_N(o.reference_N), reference_r(0)
{
  memcpy(box_length,o.box_length,3*sizeof(double));
  memcpy(box_angles,o.box_angles,3*sizeof(double));
  if (o.id) {id=new short[N]; memcpy(id,o.id,N*sizeof(short));}
  if (o.type) {type=new short[N]; memcpy(type,o.type,N*sizeof(short));}
  if (o.flags) {flags=new short[N]; memcpy(flags,o.flags,N*sizeof(short));}
  if (o.r) {r=new double[N][3]; memcpy(r,o.r,3*N*sizeof(double));}
  if (o.v) {v=new double[N][3]; memcpy(v,o.v,3*N*sizeof(double));}
  if (o.a) {a=new double[N][3]; memcpy(a,o.a,3*N*sizeof(double));}
  if (o.reference_r) {reference_r=new double[N][3];
    memcpy(reference_r,o.reference_r,3*N*sizeof(double));}
}

OLconfiguration& OLconfiguration::operator=(const OLconfiguration& c)
{
  if (this==&c) return *this;
  bool Nchanged= N!=c.N;
  N=c.N;
  name=c.name;
  time=c.time;
  step=c.step;
  memcpy(box_length,c.box_length,3*sizeof(double));
  memcpy(box_angles,c.box_angles,3*sizeof(double));
  reference_N=c.reference_N;

  if (c.id==0 || Nchanged) {delete[] id; id=0;}
  if (c.id) {
    if (id==0) id=new short[N];
    memcpy(id,c.id,N*sizeof(short));
  }
  if (c.type==0 || Nchanged) {delete[] type; type=0;}
  if (c.type) {
    if (type==0) type=new short[N];
    memcpy(type,c.type,N*sizeof(short));
  }
  if (c.flags==0 || Nchanged) {delete[] flags; flags=0;}
  if (c.flags) {
    if (flags==0) flags=new short[N];
    memcpy(flags,c.flags,N*sizeof(short));
  }
  if (c.r==0 || Nchanged) {delete[] r; r=0;}
  if (c.r) {
    if (r==0) r=new double[N][3];
    memcpy(r,c.r,3*N*sizeof(double));
  }
  if (c.v==0 || Nchanged) {delete[] v; v=0;}
  if (c.v) {
    if (v==0) v=new double[N][3];
    memcpy(v,c.v,3*N*sizeof(double));
  }
  if (c.a==0 || Nchanged) {delete[] a; a=0;}
  if (c.a) {
    if (a==0) a=new double[N][3];
    memcpy(a,c.a,3*N*sizeof(double));
  }
  if (c.reference_r==0 || Nchanged) {delete[] reference_r; reference_r=0;}
  if (c.reference_r) {
    if (reference_r==0) reference_r=new double[N][3];
    memcpy(reference_r,c.reference_r,3*N*sizeof(double));
  }
  return *this;
}

OLconfiguration& OLconfiguration::warmcpy(const OLconfiguration& c)
{
  if (this==&c) return *this;
  name=c.name;
  time=c.time;
  step=c.step;
  memcpy(box_length,c.box_length,3*sizeof(double));
  memcpy(box_angles,c.box_angles,3*sizeof(double));

  if (id && c.id) memcpy(id,c.id,N*sizeof(short));
  if (type && c.type) memcpy(type,c.type,N*sizeof(short));
  if (flags && c.flags) memcpy(flags,c.flags,N*sizeof(short));
  if (r && c.r) memcpy(r,c.r,3*N*sizeof(double));
  if (v && c.v) memcpy(v,c.v,3*N*sizeof(double));
  if (a && c.a) memcpy(a,c.a,3*N*sizeof(double));
  return *this;
}

OLconfiguration& OLconfiguration::swap(OLconfiguration& c)
{
  std::swap(N,c.N);
  std::swap(name,c.name);
  std::swap(time,c.time);
  std::swap(step,c.step);
  std::swap<double,3>(box_length,c.box_length);
  std::swap<double,3>(box_angles,c.box_angles);
  std::swap(id,c.id);
  std::swap(type,c.type);
  std::swap(flags,c.flags);
  std::swap(r,c.r);
  std::swap(v,c.v);
  std::swap(a,c.a);
}

/*
 * Init, load, and save
 *
 */
void OLconfiguration::init(const char* s)
{
  if (s==0) throw glsim::Unimplemented("Create config from scratch",HERE);
  load(s);
  step=0;
  time=0;
}

/**
This reads from a file than can contain one or many configurations, in
the latter case reads only the first.
*/
void OLconfiguration::load(const char* s) {
  OLconfig_file f(this);
  f.open(s);
  f.read_header();
  f.read_record(0);
}

void OLconfiguration::save(const char* s)
{
  OLconfig_file f(this);  // Default argument makes all fields header
  f.create(s);
  f.write_header();
}

/*
 * Periodic boundary conditions
 *
 */

void OLconfiguration::fold_coordinates()
{
  for (int n=0; n<N; n++) {
    r[n][0]-=box_length[0]*floor(r[n][0]/box_length[0]);
    r[n][1]-=box_length[1]*floor(r[n][1]/box_length[1]);
    r[n][2]-=box_length[2]*floor(r[n][2]/box_length[2]);
  }
}

/*inline*/ void OLconfiguration::fold_one(double x[])
{
  x[0]-=box_length[0]*floor(x[0]/box_length[0]);
  x[1]-=box_length[1]*floor(x[1]/box_length[1]);
  x[2]-=box_length[2]*floor(x[2]/box_length[2]);
}

/** 
  Unfold all coordinates using an internal reference.  On the first
  call, positions are not changed but a private copy of them is made
  to use as reference on the next call.  On successive calls, the
  current positions are unfolded (calling unfold_coordinates(r)) wrt
  the private reference, and the private reference is updated to the
  new positions.  In this way, an unfolded trajectory can be
  reconstructed by successively reading and unfolding.  Of course this
  will fail if the trajectory contains snapshots recorded too far
  apart, or if the simulation involves nonlocal moves.

  Call unfold_coordinates_reset() to reset the reference configuration
  (e.g. upon loading a completely different configuration).

  The stored private reference is copied upon copy construction, but
  invalidated on assingment.
 */
void OLconfiguration::unfold_coordinates()
{
  if (reference_r!=0 && reference_N==N) {
    // We have a valid private reference, unfold and updata
    unfold_coordinates(reference_r);
    memcpy(reference_r,r,3*N*sizeof(double));
  } else {
    // No valid private reference, initialize it
    unfold_coordinates_reset();
  }
}

/**
  Set current configuration as reference for unfolding, discarding the
  previous stored reference.  Will be called automatically by
  unfold_coordinates() if no reference is found.
*/
void OLconfiguration::unfold_coordinates_reset()
{
  if (reference_r==0 || reference_N!=N) {
    // Initialize reference array
    reference_N=N;
    delete[] reference_r;
    reference_r=new double[N][3];
  }

  // Update private reference
  memcpy(reference_r,r,3*N*sizeof(double));
}

/**
 Call this function to unfold all coordinates (i.e. undo previous
 folding into the periodic box done with fold_coordinates()).  This
 works by assuming the new configuration has been obtained by
 displacing the given reference positions by a small amount (e.g. by
 doing a few molecular dynamic steps).  It won't work if the reference
 positions are too far apart (in simulation time) from the current
 configuration, or if the simulation does nonlocal moves.

 \param[in]  ref  Pointer to reference positions.

 */
void OLconfiguration::unfold_coordinates(double (*ref)[3])
{
  for (int i=0; i<N; i++) {
    r[i][0] = ref[i][0] + ddiff(r[i][0],ref[i][0],box_length[0]);
    r[i][1] = ref[i][1] + ddiff(r[i][1],ref[i][1],box_length[1]);
    r[i][2] = ref[i][2] + ddiff(r[i][2],ref[i][2],box_length[2]);
  }
}

/*
 * Configuration properties
 *
 */
size_t OLconfiguration::NTypes() const
{
  std::unordered_map<short,int> types;
  for (int i=0; i<N; ++i)
    types[type[i]]++;
  return types.size();
}

int OLconfiguration::count_type(short T) const
{
  if (!type) return -1;
  int c=0;
  for (int i=0; i<N; ++i)
    if (type[i]==T) c++;
  return c;
}


/**
Returns a vector of doubles with the position of the center of mass.
If the configuration defines a type for each particle (through the
type array), you can give the masses of each type, otherwise all are
assumed to be identical.

\param[in] masses Pointer to an array of masses for each type, so that
                  `masses[type[i]]` gives the mass of particle `i`.
                  If equal to 0, particles are assumed of equal mass.
*/
std::vector<double> OLconfiguration::center_of_mass(double *masses)
{
  std::vector<double> CM;
  CM.resize(3,0.);

  if (masses) {
    double mass=0;
    for (int i=0; i<N; ++i) {
      mass+=masses[type[i]];
      for (int j=0; j<3; ++j) CM[j]+=masses[type[i]]*r[i][j];
    }
    for (int j=0; j<3; ++j) CM[j]/=mass;
  } else {
    for (int i=0; i<N; ++i)
      for (int j=0; j<3; ++j) CM[j]+=r[i][j];
  }

  return CM;
}

/*
 * cubic symmetries
 *
 */

const int OLconfiguration::CubicMatrices[48][3][3]={
  { { -1 , 0 , 0 }, { 0 , -1 , 0 }, { 0 , 0 , -1 } },
  { { -1 , 0 , 0 }, { 0 , -1 , 0 }, { 0 , 0 , 1 } },
  { { -1 , 0 , 0 }, { 0 , 0 , -1 }, { 0 , -1 , 0 } },
  { { -1 , 0 , 0 }, { 0 , 0 , 1 },  { 0 , -1 , 0 } },
  { { -1 , 0 , 0 }, { 0 , 0 , -1 }, { 0 , 1 , 0 } },
  { { -1 , 0 , 0 }, { 0 , 0 , 1 }, { 0 , 1 , 0 } },
  { { -1 , 0 , 0 }, { 0 , 1 , 0 }, { 0 , 0 , -1 } },
  { { -1 , 0 , 0 }, { 0 , 1 , 0 }, { 0 , 0 , 1 } },
  { { 0 , -1 , 0 }, { -1 , 0 , 0 }, { 0 , 0 , -1 } },
  { { 0 , -1 , 0 }, { -1 , 0 , 0 }, { 0 , 0 , 1 } },
  { { 0 , 0 , -1 }, { -1 , 0 , 0 }, { 0 , -1 , 0 } },
  { { 0 , 0 , 1 }, { -1 , 0 , 0 }, { 0 , -1 , 0 } },
  { { 0 , 0 , -1 }, { -1 , 0 , 0 }, { 0 , 1 , 0 } },
  { { 0 , 0 , 1 }, { -1 , 0 , 0 }, { 0 , 1 , 0 } },
  { { 0 , 1 , 0 }, { -1 , 0 , 0 }, { 0 , 0 , -1 } }, 
  { { 0 , 1 , 0 }, { -1 , 0 , 0 },{ 0 , 0 , 1 } },
  { { 0 , -1 , 0 }, { 0 , 0 , -1 },{ -1 , 0 , 0 } },
  { { 0 , -1 , 0 }, { 0 , 0 , 1 }, { -1 , 0 , 0 } },
  { { 0 , 0 , -1 }, { 0 , -1 , 0 }, { -1 , 0 , 0 } },
  { { 0 , 0 , 1 }, { 0 , -1 , 0 }, { -1 , 0 , 0 } },
  { { 0 , 0 , -1 }, { 0 , 1 , 0 }, { -1 , 0 , 0 } },
  { { 0 , 0 , 1 }, { 0 , 1 , 0 }, { -1 , 0 , 0 } },
  { { 0 , 1 , 0 }, { 0 , 0 , -1 }, { -1 , 0 , 0 } },
  { { 0 , 1 , 0 }, { 0 , 0 , 1 }, { -1 , 0 , 0 } },
  { { 0 , -1 , 0 }, { 0 , 0 , -1 }, { 1 , 0 , 0 } },
  { { 0 , -1 , 0 }, { 0 , 0 , 1 }, { 1 , 0 , 0 } },
  { { 0 , 0 , -1 }, { 0 , -1 , 0 }, { 1 , 0 , 0 } },
  { { 0 , 0 , 1 }, { 0 , -1 , 0 }, { 1 , 0 , 0 } },
  { { 0 , 0 , -1 }, { 0 , 1 , 0 }, { 1 , 0 , 0 } },
  { { 0 , 0 , 1 }, { 0 , 1 , 0 }, { 1 , 0 , 0 } },
  { { 0 , 1 , 0 }, { 0 , 0 , -1 }, { 1 , 0 , 0 } },
  { { 0 , 1 , 0 }, { 0 , 0 , 1 }, { 1 , 0 , 0 } },
  { { 0 , -1 , 0 }, { 1 , 0 , 0 }, { 0 , 0 , -1 } },
  { { 0 , -1 , 0 }, { 1 , 0 , 0 }, { 0 , 0 , 1 } },
  { { 0 , 0 , -1 }, { 1 , 0 , 0 }, { 0 , -1 , 0 } },
  { { 0 , 0 , 1 }, { 1 , 0 , 0 }, { 0 , -1 , 0 } },
  { { 0 , 0 , -1 }, { 1 , 0 , 0 }, { 0 , 1 , 0 } },
  { { 0 , 0 , 1 }, { 1 , 0 , 0 }, { 0 , 1 , 0 } },
  { { 0 , 1 , 0 }, { 1 , 0 , 0 }, { 0 , 0 , -1 } },
  { { 0 , 1 , 0 }, { 1 , 0 , 0 }, { 0 , 0 , 1 } },
  { { 1 , 0 , 0 }, { 0 , -1 , 0 }, { 0 , 0 , -1 } },
  { { 1 , 0 , 0 }, { 0 , -1 , 0 }, { 0 , 0 , 1 } },
  { { 1 , 0 , 0 }, { 0 , 0 , -1 }, { 0 , -1 , 0 } },
  { { 1 , 0 , 0 }, { 0 , 0 , 1 }, { 0 , -1 , 0 } },
  { { 1 , 0 , 0 }, { 0 , 0 , -1 }, { 0 , 1 , 0 } },
  { { 1 , 0 , 0 }, { 0 , 0 , 1 }, { 0 , 1 , 0 } },
  { { 1 , 0 , 0 }, { 0 , 1 , 0 }, { 0 , 0 , -1 } },
  { { 1 , 0 , 0 }, { 0 , 1 , 0 }, { 0 , 0 , 1 } }
  } ;

void apply_cubic_operation(OLconfiguration& conf,
			   double (*r)[3],int n)
{
  const int (*m)[3]=conf.CubicMatrices[n];
  double out[3];
  
  for (int i=0; i<conf.N; i++) {
    out[0]=m[0][0]*conf.r[i][0] + m[0][1]*conf.r[i][1] + m[0][2]*conf.r[i][2];
    out[1]=m[1][0]*conf.r[i][0] + m[1][1]*conf.r[i][1] + m[1][2]*conf.r[i][2];
    out[2]=m[2][0]*conf.r[i][0] + m[2][1]*conf.r[i][1] + m[2][2]*conf.r[i][2];
    memcpy(r[i],out,3*sizeof(double));
  }
}



/*****************************************************************************
 *
 * OLconfig_file
 *
 */

OLconfig_file::OLconfig_file(OLconfiguration *buffer,
			     OLconfig_file::options o) :
  HDF_record_file(1),
  version(1),
  cbuffer(buffer),
  opt(o)
{
  if (cbuffer==0) {
    cbuffer=new OLconfiguration();
    own_buffer=true;
  } else
    own_buffer=false;
}

OLconfig_file::~OLconfig_file()
{
  if (own_buffer) delete cbuffer;
}

void OLconfig_file::create(const char* f)
{
  fname=f;
  HDF_record_file::open(fname.c_str(),f_replace,"OLconfig_file",version_exact,
			version);
}

void OLconfig_file::open(const char *f)
{
  fname=f;
  HDF_record_file::open(fname.c_str(),f_append,"OLconfig_file",version_exact,
			version);
}

void OLconfig_file::open_ro(const char *f)
{
  fname=f;
  HDF_record_file::open(fname.c_str(),f_readonly,"OLconfig_file",version_exact,
			version);
}

void OLconfig_file::declare_header_fields(mode fmode)
{
  declare_attribute("N",&(cbuffer->N));
  declare_field(f_header,"name",&(cbuffer->name));
  if (!opt.time_frame_) {
    declare_field(f_header,"time",&(cbuffer->time));
    declare_field(f_header,"step",&(cbuffer->step));
  }
  if (!opt.box_frame_) {
    declare_field(f_header,"box_length",cbuffer->box_length,3);
    declare_field(f_header,"box_angles",cbuffer->box_angles,3);
  }

  if (fmode==f_replace) {
    if (!opt.id_frame_)
      create_if_present("id","has_id",cbuffer->id,f_header);
    if (!opt.type_frame_)
      create_if_present("type","has_type",cbuffer->type,f_header);
    if (!opt.flags_frame_)
      create_if_present("flags","has_flags",cbuffer->flags,f_header);
    if (!opt.r_frame_)
      create_if_present("r","has_r",cbuffer->r,f_header);
    if (!opt.v_frame_)
      create_if_present("v","has_v",cbuffer->v,f_header);
    if (!opt.a_frame_)
      create_if_present("a","has_a",cbuffer->a,f_header);
  } else {
    if (!opt.id_frame_)
      open_if_present("id","has_id",cbuffer->id);
    if (!opt.id_frame_)
      open_if_present("type","has_type",cbuffer->type);
    if (!opt.id_frame_)
      open_if_present("flags","has_flags",cbuffer->flags);
    if (!opt.id_frame_)
      open_if_present("r","has_r",cbuffer->r);
    if (!opt.id_frame_)
      open_if_present("v","has_v",cbuffer->v);
    if (!opt.id_frame_)
      open_if_present("a","has_a",cbuffer->a);
  }
}

void OLconfig_file::declare_record_fields(mode fmode)
{
  if (opt.time_frame_) {
    declare_field(f_record,"time",&(cbuffer->time));
    declare_field(f_record,"step",&(cbuffer->step));
  }
  if (opt.box_frame_) {
    declare_field(f_record,"box_length",cbuffer->box_length,3);
    declare_field(f_record,"box_angles",cbuffer->box_angles,3);
  }
  if (fmode==f_replace) {
    if (opt.id_frame_)
      create_if_present("id","has_id",cbuffer->id,f_record);
    if (opt.type_frame_)
      create_if_present("type","has_type",cbuffer->type,f_record);
    if (opt.flags_frame_)
      create_if_present("flags","has_flags",cbuffer->flags,f_record);
    if (opt.r_frame_)
      create_if_present("r","has_r",cbuffer->r,f_record);
    if (opt.v_frame_)
      create_if_present("v","has_v",cbuffer->v,f_record);
    if (opt.a_frame_)
      create_if_present("a","has_a",cbuffer->a,f_record);
  }
}


} /* namespace */
