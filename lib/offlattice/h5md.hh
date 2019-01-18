/*
 * h5md.hh -- class interface for H5MD trajectory files
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

#ifndef H5MD_HH
#define H5MD_HH

#include "hdf_file.hh"
#include "olconfiguration.hh"

namespace glsim {

/** \class H5MD
    \ingroup Offlatticeconf

This is a helper class to write `OLconfiguration`screate trajectory
files in H5MD format.  So far it is rather primitive, writing only
positions vs time and with a fixed box.
*/
class H5MD : public H5file {
public:
  H5MD() :
    position_step_t(1),
    position_time_t(1),
    position_value_t(1),
    last_record(-1)
  {}
  void create(const char *fname,int Nparticles,const char *author,
	      const char *author_email=0);
  void set_box(double b[]);
  void append_positions(const OLconfiguration &c);

private:
  H5group h5md,particles;
  H5group author,creator;
  H5group all,box,position;

  H5set                   position_step,position_time,position_value;
  H5simple_type<long,1>   position_step_t;
  H5simple_type<double,1> position_time_t,position_value_t;

  hsize_t                 last_record;
} ;

} /* namespace */


#endif /* H5MD_HH */
