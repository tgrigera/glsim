/*
 * h5md.cc -- class interface for H5MD trajectory files (definitions)
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

#include "../config.h"

#include "h5md.hh"

namespace glsim {

void H5MD::create(const char *fname,int Nparticles,const char *author_name,
		  const char *author_email)
{
  H5file::create(fname,H5F_ACC_TRUNC);

  h5md.create(this->id(),"h5md");
  int v[2]={1,0};
  HE(H5LTset_attribute_int(h5md.id(),".","version",v,2));

  author.create(h5md.id(),"author");
  HE(H5LTset_attribute_string(author.id(),".","name",author_name));
  if (author_email!=0)
    HE(H5LTset_attribute_string(author.id(),".","email",author_email));

  creator.create(h5md.id(),"creator");
  HE(H5LTset_attribute_string(creator.id(),".","name","glsim"));
  HE(H5LTset_attribute_string(creator.id(),".","version",PACKAGE_VERSION));

  particles.create(this->id(),"particles");
  all.create(particles.id(),"all");

  box.create(all.id(),"box");
  int d=3;
  hid_t H5_maj=H5Acreate2(box.id(),"dimension",
			  H5T_NATIVE_INT,H5_scalar_s,H5P_DEFAULT,H5P_DEFAULT);
  HE(H5_maj);
  HE(H5Awrite(H5_maj,H5T_NATIVE_INT,&d));
  HE(H5Aclose(H5_maj));
  char boundary[3][9]={"periodic","periodic","periodic"};
  H5simple_type<char,1> ST(9);
  hsize_t cur[1]={3};
  hsize_t max[1]={3};
  hid_t H5_arr3=H5Screate_simple(1,cur,max);
  HE(H5_arr3);
  H5_maj=H5Acreate2(box.id(),"boundary",ST.id(),
		    H5_arr3,H5P_DEFAULT,H5P_DEFAULT);
  HE(H5_maj);
  HE(H5Awrite(H5_maj,ST.id(),boundary[0]));
  HE(H5Aclose(H5_maj));
  HE(H5Sclose(H5_arr3));

  position.create(all.id(),"position");

  // datasets in position group
  hsize_t chunks[3]={1,0,0};
  hid_t PL=H5Pcreate(H5P_DATASET_CREATE);
  HE(PL);
  HE(H5Pset_chunk(PL,1,chunks));
  cur[0]=0;
  max[0]=H5S_UNLIMITED;
  hid_t H5_frame_s;
  HE(H5_frame_s=H5Screate_simple(1,cur,max));

  position_step.create(position.id(),"step",position_step_t.id(),H5_frame_s,
		       H5P_DEFAULT,PL);
  position_time.create(position.id(),"time",position_time_t.id(),H5_frame_s,
		       H5P_DEFAULT,PL);

  H5space_simple<3> pds;
  pds.create(0,Nparticles,3,H5S_UNLIMITED,Nparticles,3);
  chunks[0]=1;
  chunks[1]=Nparticles;
  chunks[2]=3;
  HE(H5Pset_chunk(PL,3,chunks));
  position_value.create(position.id(),"value",position_value_t.id(),pds.id(),
			H5P_DEFAULT,PL);
  HE(H5Sclose(H5_frame_s));
  HE(H5Pclose(PL));
}

void H5MD::set_box(double b[])
{
  H5simple_type<double,1> boxt(1);
  H5space_simple<1> aspace;
  aspace.create(3,3);
  H5set edges;
  edges.create(box.id(),"edges",boxt.id(),aspace.id());
  HE(H5Dwrite(edges.id(),boxt.id(),H5S_ALL,H5S_ALL,H5P_DEFAULT,
    		b));
}

void H5MD::append_positions(const OLconfiguration &conf)
{
  // Define the slab (just one point) for i/o
  last_record++;
  hsize_t start[3] = {last_record,0,0};
  hsize_t count[3] = {1,(unsigned) conf.N,3};

  // resize
  hsize_t newsize[3];
  newsize[0]=last_record+1;
  newsize[1]=conf.N;
  newsize[2]=3;
  HE(H5Dset_extent(position_step.id(),newsize));
  HE(H5Dset_extent(position_time.id(),newsize));
  HE(H5Dset_extent(position_value.id(),newsize));

  // First need to select the record through a hyperslab
  // (need a copy of the dataspace)
  hid_t selection;
  HE(selection=H5Dget_space(position_step.id()));
  HE(H5Sselect_hyperslab(selection,H5S_SELECT_SET,
			 start,0,count,0));   
  HE(H5Dwrite(position_step.id(),position_step_t.id(),
	      H5_scalar_s,selection,H5P_DEFAULT,&conf.step));
  HE(H5Sclose(selection));
  HE(selection=H5Dget_space(position_time.id()));
  HE(H5Sselect_hyperslab(selection,H5S_SELECT_SET,
			 start,0,count,0));   
  HE(H5Dwrite(position_time.id(),position_time_t.id(),
	      H5_scalar_s,selection,H5P_DEFAULT,&conf.time));
  HE(H5Sclose(selection));
  HE(selection=H5Dget_space(position_value.id()));
  HE(H5Sselect_hyperslab(selection,H5S_SELECT_SET,
			 start,0,count,0));   
  H5space_simple<2> wpos;
  wpos.create(conf.N,3,conf.N,3);
  HE(H5Dwrite(position_value.id(),position_value_t.id(),
	      wpos.id(),selection,H5P_DEFAULT,conf.r));
  HE(H5Sclose(selection));
}



} /* namespace */
