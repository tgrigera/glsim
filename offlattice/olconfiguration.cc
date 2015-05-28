/*
 * olconfiguration.hh -- definitions for class OLconfiguration
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

#include <string.h>
#include <math.h>
#include <utility>

#include "cerrors.h"
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
  N=0;
}

OLconfiguration::OLconfiguration(const OLconfiguration& o) :
  Configuration(o.name),
  N(o.N),
  time(o.time),
  step(o.step),
  id(0), type(0), flags(0), r(0), v(0), a(0)
{
  memcpy(box_length,o.box_length,3*sizeof(double));
  memcpy(box_angles,o.box_angles,3*sizeof(double));
  if (o.id) {id=new short[N]; memcpy(id,o.id,N*sizeof(short));}
  if (o.type) {type=new short[N]; memcpy(type,o.type,N*sizeof(short));}
  if (o.flags) {flags=new short[N]; memcpy(flags,o.flags,N*sizeof(short));}
  if (o.r) {r=new double[N][3]; memcpy(r,o.r,3*N*sizeof(double));}
  if (o.v) {v=new double[N][3]; memcpy(v,o.v,3*N*sizeof(double));}
  if (o.a) {a=new double[N][3]; memcpy(a,o.a,3*N*sizeof(double));}
}

OLconfiguration& OLconfiguration::operator=(const OLconfiguration& c)
{
  if (this==&c) return *this;
  N=c.N;
  name=c.name;
  time=c.time;
  step=c.step;
  memcpy(box_length,c.box_length,3*sizeof(double));
  memcpy(box_angles,c.box_angles,3*sizeof(double));
  delete[] id; id=0;
  delete[] type; type=0;
  delete[] flags; flags=0;
  delete[] r; r=0;
  delete[] v; v=0;
  delete[] a; a=0;
  if (c.id) {id=new short[N]; memcpy(id,c.id,N*sizeof(short));}
  if (c.type) {type=new short[N]; memcpy(type,c.type,N*sizeof(short));}
  if (c.flags) {flags=new short[N]; memcpy(flags,c.flags,N*sizeof(short));}
  if (c.r) {r=new double[N][3]; memcpy(r,c.r,3*N*sizeof(double));}
  if (c.v) {v=new double[N][3]; memcpy(v,c.v,3*N*sizeof(double));}
  if (c.a) {a=new double[N][3]; memcpy(a,c.a,3*N*sizeof(double));}
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
}

/**
This reads from a file than can contain one or many configurations, in
the latter case reads only the first.
*/
void OLconfiguration::load(const char* s) {
  OLconfig_file f(s,this);
  f.open();
  f.read_header();
  f.read_record(0);
}

void OLconfiguration::save(const char* s)
{
  OLconfig_file f(s,this);  // Default argument makes all fields header
  f.create();
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

/*****************************************************************************
 *
 * OLconfig_file
 *
 */

OLconfig_file::OLconfig_file(const char *f,OLconfiguration *buffer) :
  HDF_record_file(1),
  version(1),
  fname(f),
  cbuffer(buffer)
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

void OLconfig_file::create()
{
  HDF_record_file::open(fname.c_str(),f_replace,"OLconfig_file",version_exact,
			version);
}

void OLconfig_file::open()
{
  HDF_record_file::open(fname.c_str(),f_append,"OLconfig_file",version_exact,
			version);
}

void OLconfig_file::declare_header_fields(mode fmode)
{
  declare_attribute("N",&(cbuffer->N));
  declare_field(f_header,"name",&(cbuffer->name));
  declare_field(f_header,"time",&(cbuffer->time));
  declare_field(f_header,"step",&(cbuffer->step));
  declare_field(f_header,"box_length",cbuffer->box_length,3);
  declare_field(f_header,"box_angles",cbuffer->box_angles,3);

  if (fmode==f_replace) {
    create_if_present("id","has_id",cbuffer->id);
    create_if_present("type","has_type",cbuffer->type);
    create_if_present("flags","has_flags",cbuffer->flags);
    create_if_present("r","has_r",cbuffer->r);
    create_if_present("v","has_v",cbuffer->v);
    create_if_present("a","has_a",cbuffer->a);
  } else {
    open_if_present("id","has_id",cbuffer->id);
    open_if_present("type","has_type",cbuffer->type);
    open_if_present("flags","has_flags",cbuffer->flags);
    open_if_present("r","has_r",cbuffer->r);
    open_if_present("v","has_v",cbuffer->v);
    open_if_present("a","has_a",cbuffer->a);
  }
}

void OLconfig_file::declare_record_fields(mode fmode)
{}


} /* namespace */
