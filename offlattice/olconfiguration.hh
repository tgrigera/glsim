/*
 * olconfiguration.hh -- declaration of class OLconfiguration
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

#ifndef OLCONFIGURATION_HH
#define OLCONFIGURATION_HH

#include "configuration.hh"
#include "exception.hh"
#include "hdf_file.hh"

namespace glsim {

/*****************************************************************************/

/** \class OLconfiguration
    \ingroup Offlattice

*/
class OLconfiguration : public Configuration {
public:
  int       N;
  double    time;
  long      step;
  double    box_length[3];
  double    box_angles[3];

  short     *id;
  short     *type;
  short     *flags;
  double    (*r)[3];
  double    (*v)[3];
  double    (*a)[3];
  
public:
  OLconfiguration();
  OLconfiguration(const std::string &title);
  OLconfiguration(const OLconfiguration&);
  OLconfiguration& operator=(const OLconfiguration&);
  ~OLconfiguration();
  OLconfiguration& swap(OLconfiguration& c);

  using Configuration::init;
  using Configuration::load;
  using Configuration::save;
  void load(const char* fname);
  void save(const char* fname);
  void init(const char* fname);

  // void fold_coordinates();
  // void fold_one(double x[]);

  // double distancesq(const double x[],const double y[]) const;
  // double distancesq(int i,int j) const;

  // int CountType(short t) const;
  // int CountFlag(short f) const;

  // void ApplyCubicSymmetry(OLconfiguration&,int n);
} ;

inline OLconfiguration::OLconfiguration() :
  Configuration(""), 
  N(0),
  time(0),
  step(0),
  id(0), type(0), flags(0), r(0), v(0), a(0)
{
  box_length[0]=box_length[1]=box_length[2]=0;
  box_angles[0]=box_angles[1]=box_angles[2]=0;
}

inline OLconfiguration::OLconfiguration(const std::string &title) :
  Configuration(title), 
  N(0),
  time(0),
  step(0),
  id(0), type(0), flags(0), r(0), v(0)
{
  box_length[0]=box_length[1]=box_length[2]=0;
  box_angles[0]=box_angles[1]=box_angles[2]=0;
}

/*****************************************************************************/

/** \class OLconfig_file
    \ingroup Offlattice
*/
class OLconfig_file : public HDF_record_file {
public:
  OLconfig_file(const char *fname,OLconfiguration *buffer=0);
  ~OLconfig_file();
  
  void create();
  void open();

private:
  const int       version;
  std::string     fname;
  OLconfiguration *cbuffer;
  bool            own_buffer;

  void declare_header_fields(mode);
  void declare_record_fields(mode);
  template <typename FTYPE>
  void create_if_present(const char* field_name,const char* aname,FTYPE* p);
  template <typename FTYPE>
  void open_if_present(const char* field_name,const char* aname,FTYPE* (&p));
} ;

template <typename FTYPE>
void OLconfig_file::create_if_present(const char* field_name,const char* aname,
				      FTYPE* p)
{
  bool hasp;
  hasp= p!=0;
  declare_attribute(aname,&hasp);
  if (hasp)
    declare_field(f_header,field_name,p,cbuffer->N);
}

template <typename FTYPE>
void OLconfig_file::open_if_present(const char* field_name,const char* aname,
				    FTYPE* (&p))
{
  bool hasp;
  declare_attribute(aname,&hasp);
  if (hasp) {
    delete[] p;
    p=new FTYPE[cbuffer->N];
    declare_field(f_header,field_name,p,cbuffer->N);
  } else {
    delete[] p;
    p=0;
  }
}

} /* namespace */

#endif /* OLCONFIGURATION_HH */
