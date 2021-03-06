/*
 * olconfiguration.hh -- declaration of class OLconfiguration
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

#ifndef OLCONFIGURATION_HH
#define OLCONFIGURATION_HH

#include <cmath>

#include "configuration.hh"
#include "exception.hh"
#include "hdf_file.hh"

namespace glsim {

/*****************************************************************************/

/** \class OLconfiguration
    \ingroup OfflatticeCONF
    \brief A configuration class for off-lattice systems

*/
class OLconfiguration : public Configuration {
public:
  int       N;             ///< Number of particles
  double    time;          ///< Time at which the configuration was recorded
  long      step;          ///< Step number at which the configuration was recorded
  double    box_length[3]; ///< Length of the box vectors
  double    box_angles[3]; ///< Angles

  short     *id;           ///< Particle (unique) identification
  short     *type;         ///< Particle type (use as index to find properties)
  short     *flags;        ///< Flags for special signalling
  double    (*r)[3];       ///< Positions
  double    (*v)[3];       ///< Velocities
  double    (*a)[3];       ///< Accelerations or forces
  
public:
  /// \name Creation, destruction, and copying
  /// @{
  OLconfiguration();
  OLconfiguration(const std::string &title);
  OLconfiguration(const OLconfiguration&);
  /// Copy the whole object, including the reference configuration
  /// *Warning:* This may change the positions of the buffers (arrays) if the source
  /// configuration has a different N or contains a different set of nonempty fields
  OLconfiguration& operator=(const OLconfiguration&);
  /// Copy all the data except reference configuration, assuming that all buffers
  /// are correctly allocated (same size, same set of nonempty fields).  Undefined behavior
  /// will result if N or nonempty fields differ, or if buffers were not previously allocate
  OLconfiguration& warmcpy(const OLconfiguration&);
  ~OLconfiguration();
  OLconfiguration& swap(OLconfiguration& c);

  /// @}
  /// \name Initialization, load and save
  /// @{
  using Configuration::init;
  using Configuration::load;
  using Configuration::save;
  void load(const char* fname);
  void save(const char* fname);
  void init(const char* fname);

  /// @}
  /// \name Periodic boundary conditions
  /// @{
  void fold_coordinates();     ///< Move all particles inside the primary box (shifting by a multiple of the box length)
  void unfold_coordinates(double (*ref)[3]); ///< Undo periodic folding, using a proided reference configuration
  void unfold_coordinates(); ///< Undo periodic folding, using internal reference configuration
  void unfold_coordinates_reset(); ///< Reset reference configuration for unfold_coordinates()
  void fold_one(double x[]);   ///< Apply PBCs to one particle

  double ddiff(const double a,const double b,const double box_length) const;
  double distancesq(const double x[],const double y[]) const;
  double distancesq(int i,int j) const;

  /// @}
  /// \name Configuration properties
  /// @{
  size_t NTypes() const;
  int    count_type(short t) const;
  // int CountFlag(short f) const;
  double volume() const;  ///< Box volume
  double number_density() const;  ///< conf.N / volume
  std::vector<double> center_of_mass(double *masses=0);  ///< Compute the center of mass

  /// @}

private:
  int       reference_N;
  double    (*reference_r)[3];

  static const int CubicMatrices[48][3][3];
  friend void apply_cubic_operation(OLconfiguration& conf,double (*r)[3],int n);

} ;

/** \brief  Apply cubic symmetry operation

Apply the nth (from 0 to 47) cubic symmetry operation to given
configuration.  Output coordinates are put in r (may coincide with
conf.r)
*/
void apply_cubic_operation(OLconfiguration& conf,double (*r)[3],int n);

inline OLconfiguration::OLconfiguration() :
  Configuration(""), 
  N(0),
  time(0),
  step(0),
  id(0), type(0), flags(0), r(0), v(0), a(0),
  reference_N(0), reference_r(0)
{
  box_length[0]=box_length[1]=box_length[2]=0;
  box_angles[0]=box_angles[1]=box_angles[2]=0;
}

inline OLconfiguration::OLconfiguration(const std::string &title) :
  Configuration(title), 
  N(0),
  time(0),
  step(0),
  id(0), type(0), flags(0), r(0), v(0),
  reference_N(0), reference_r(0)
{
  box_length[0]=box_length[1]=box_length[2]=0;
  box_angles[0]=box_angles[1]=box_angles[2]=0;
}

/* 
 * Functions for periodic boundary conditions
 *
 */

/** Returns a-b shifted by periodic images as necessary so that |a-b|<box/2
 */
inline double OLconfiguration::ddiff(const double a,const double b,
				     const double box_length) const
{
  double temp=a-b;
  return temp-box_length*round(temp/box_length);
}

inline double OLconfiguration::distancesq(const double x[],
					  const double y[]) const
{
  double dx=ddiff(x[0],y[0],box_length[0]);
  double dy=ddiff(x[1],y[1],box_length[1]);
  double dz=ddiff(x[2],y[2],box_length[2]);
  return dx*dx+dy*dy+dz*dz;
}

inline double OLconfiguration::distancesq(int i,int j) const
{
  double dx=ddiff(r[i][0],r[j][0],box_length[0]);
  double dy=ddiff(r[i][1],r[j][1],box_length[1]);
  double dz=ddiff(r[i][2],r[j][2],box_length[2]);
  return dx*dx+dy*dy+dz*dz;
}

/*
 * Configuration properties
 *
 */

inline double OLconfiguration::number_density() const
{
  return N/volume();
}

inline double OLconfiguration::volume() const
{
  return box_length[0]*box_length[1]*box_length[2];
}


/*****************************************************************************/

/** \class OLconfig_file
    \ingroup OfflatticeCONF
*/
class OLconfig_file : public HDF_record_file {
public:
  struct options;
  OLconfig_file(OLconfiguration *buffer=0,
		OLconfig_file::options o=OLconfig_file::options());
  ~OLconfig_file();
  
  void create(const char* fname);
  void open(const char* fname);
  void open_ro(const char* fname);

  struct options {
    bool time_frame_;
    bool box_frame_;
    bool id_frame_;
    bool type_frame_;
    bool flags_frame_;
    bool r_frame_;
    bool v_frame_;
    bool a_frame_;

    options() :
      time_frame_(false), box_frame_(false), id_frame_(false),
      type_frame_(false), flags_frame_(false),
      r_frame_(false), v_frame_(false), a_frame_(false) {}

    options& time_frame() {time_frame_=true; return *this;}
    options& box_frame() {time_frame_=true; box_frame_=true; return *this;}
    options& id_frame() {time_frame_=true; id_frame_=true; return *this;}
    options& type_frame() {time_frame_=true; type_frame_=true; return *this;}
    options& flags_frame() {time_frame_=true; flags_frame_=true; return *this;}
    options& r_frame() {time_frame_=true; r_frame_=true; return *this;}
    options& v_frame() {time_frame_=true; v_frame_=true; return *this;}
    options& a_frame() {time_frame_=true; a_frame_=true; return *this;}
  } ;

private:
  const int       version;
  OLconfiguration *cbuffer;
  options         opt;
  std::string     fname;
  bool            own_buffer;

  void declare_header_fields(mode);
  void declare_record_fields(mode);
  template <typename FTYPE>
  void create_if_present(const char* field_name,const char* aname,FTYPE* p,
			 field_kind fkind);
  template <typename FTYPE>
  void open_if_present(const char* field_name,const char* aname,FTYPE* (&p));
} ;

template <typename FTYPE>
void OLconfig_file::create_if_present(const char* field_name,const char* aname,
				      FTYPE* p,field_kind fkind)
{
  bool hasp;
  hasp= p!=0;
  declare_attribute(aname,&hasp);
  if (hasp)
    declare_field(fkind,field_name,p,cbuffer->N);
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
