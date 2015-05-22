/*
 * hdf_file.cc -- class interface for HDF5 files (definitions)
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

#include "hdf_file.hh"

namespace glsim {

/******************************************************************************
 *
 * H5file
 *
 */

int   H5file::file_count=0;
hid_t H5file::H5_scalar_s;
hid_t H5file::H5_file_access_p;

H5file::H5file() :
  hfile_id(-1)
{
  if (file_count==0) {
    HE(H5open());
    HE(H5Eset_auto2(H5E_DEFAULT,0,0));
    H5_scalar_s=H5Screate(H5S_SCALAR);
    H5Pset_libver_bounds(H5_file_access_p, H5F_LIBVER_18, H5F_LIBVER_18);
  }
  file_count++;
}

H5file::~H5file()
{
  file_count--;
  if (file_count==0) {
    HE(H5Sclose(H5_scalar_s));
    HE(H5close());
  }
}

void H5file::handle_HDF_error(const Source_context c)
{
  // std::cerr << "HD5 error at " << p() << '\n';
  H5Eprint2(H5E_DEFAULT,stderr);
  throw HDF_error(c);
}


/*****************************************************************************
 *
 * HDF_record_file
 *
 */

// Do not try to overload the constructor so that it calls open directly,
// it won't work because open needs to call a virtual (init_derived) so that
// the object must be fully constructed before open is called. 
// The derived base can easily call open from the constructor 
// if it is not expected to be inherited from

/** \brief Open new or existing file

 \param  fname   File name

 \param  m       File mode, one of
                 - f_replace: create or overwrite file
		 - f_append: open existing file with read/write random access
		             similar to fopen's `a+` mode, actually allows 
			     reading or writing of header or any record
		 - f_readonly: open for reading

 \param  tit     A title for the file (stored as attribure in the root group)

 \param vr The file object and the file on disk might have different
          versions. If vr=any, disk file version is ignored, if
          vr=exact then an exception is thrown if versions don't
          match, if vr=min, exception is thrown if the disk file
          version is less thatn minver.  A disk file newer than the object 
	  always throws exception.
 
 */
void HDF_record_file::open(const char *fname,const mode& m,const char *tit,
			   version_require vr,long minver)
{
  close();
  file_mode=m;
  switch(file_mode) {
  case f_readonly:
    file=H5Fopen(fname,H5F_ACC_RDONLY,H5P_DEFAULT);
    break;
  case f_replace:
    file=H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5_file_access_p);
    break;
  case f_append:
    file=H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
    break;
  }
  HE(file);
  rootg=H5Gopen2(file,"/",H5P_DEFAULT);
  HE(rootg);
  init(tit,vr,minver);
}

void HDF_record_file::open(hid_t gr,const mode& m,const char *tit,
			   version_require vr,long minver)
{
  close();
  file_mode=m;
  switch(file_mode) {
  case f_readonly:
    break;
  case f_replace:
    throw Unimplemented("Replace when opening a datagroup",HERE);
    break;
  case f_append:
    break;
  }
  init(tit,vr,minver);
}

void HDF_record_file::init(const char* tit,version_require vr,long minver)
{
  hid_t H5_maj,H5_min;

  switch(file_mode) {
  case f_readonly:
  case f_append:
    HE(H5_maj=H5Aopen(rootg,"structure_version",H5P_DEFAULT));
    HE(H5Aread(H5_maj,H5T_NATIVE_SHORT,&disk_file_version));
    HE(H5Aclose(H5_maj));

    if (disk_file_version>file_object_version) {
      std::string m="File version (" + std::to_string(disk_file_version) +
	") is newer than code (" + std::to_string(file_object_version) +
	")";
      throw Unsupported_file_version(m,HERE);
    }
    switch (vr) {
    case version_any: break;
    case version_exact:
      if (file_object_version!=disk_file_version) {
	std::string m="Unsupported file version " +
	  std::to_string(disk_file_version) + " (required =" +
	  std::to_string(file_object_version) + ")";
	throw Unsupported_file_version(m,HERE);
      }
      break;
    case version_min:
      if (disk_file_version<minver) {
	std::string m="Unsupported file version " +
	  std::to_string(disk_file_version) + " (required >=" +
	  std::to_string(minver) +")";
	throw Unsupported_file_version(m,HERE);
      }
      break;
    }
    user_attrs.open(rootg,"User attributes");
    break;

  case f_replace:
    H5_maj=H5Acreate2(rootg,"structure_version",
    		      H5T_NATIVE_LONG,H5_scalar_s,H5P_DEFAULT,H5P_DEFAULT);
    HE(H5_maj);
    HE(H5Awrite(H5_maj,H5T_NATIVE_SHORT,&file_object_version));
    HE(H5Aclose(H5_maj));
    disk_file_version=file_object_version;
    
    // The above operations could be done with a single call to the
    // following function from the "light" interface, except that the
    // dataspace is not scalar but "simple" of size 1
    // HE(H5LTset_attribute_short(rootg,"/","structure_version_major",
    // 			       &file_version_.major,1));

    if (tit) H5LTset_attribute_string(rootg,".","title",tit);
    user_attrs.create(rootg,"User attributes");
    break;
  }
  // Create the frame (unlimited) dataspace
  hsize_t cur[1]={0};
  hsize_t max[1]={H5S_UNLIMITED};
  H5_frame_s=H5Screate_simple(1,cur,max);
  HE(H5_frame_s);

  declare_header_fields(file_mode);
  if (file_mode!=f_replace) read_header();
  declare_record_fields(file_mode);
  if (file_mode==f_replace) nrecords=0;
  else {
    if (record_fields.empty())
      nrecords=0;
    else {
      field_info fi=record_fields.front();
      hsize_t dim[1];
      HE(H5Sget_simple_extent_dims(fi.dataspace_id,dim,0));
      nrecords=dim[0];
    }
  }
}

void HDF_record_file::close()
{
  if (file<0 && rootg<0) return;
  for (const field_info &fd : header_fields) {
    if (fd.mem_datatype!=fd.file_datatype)
      delete fd.mem_datatype;
    delete fd.file_datatype;
    if (fd.dataspace_id!=H5_scalar_s)
      HE(H5Sclose(fd.dataspace_id));
    HE(H5Dclose(fd.dataset_id));
  }
  for (const field_info &fd : record_fields) {
    if (fd.mem_datatype!=fd.file_datatype)
      delete fd.mem_datatype;
    delete fd.file_datatype;
    if (fd.dataspace_id!=H5_frame_s)
      HE(H5Sclose(fd.dataspace_id));
    HE(H5Dclose(fd.dataset_id));
  }
  header_fields.clear();
  record_fields.clear();
  HE(H5Sclose(H5_frame_s));
  user_attrs.close();
  HE(H5Gclose(rootg));
  rootg=-1;
  if (file>=0) HE(H5Fclose(file));
  file=-1;
}

/*****************************************************************************
 *
 * Declare attributes
 *
 */

template <>
void HDF_record_file::declare_attribute(const char* name,bool *buffer,
					hsize_t nelem)
{
  short bval;
  
  switch (file_mode) {
  case f_readonly:
  case f_append:
    HE(H5LTget_attribute_short(user_attrs.id(),".",name,&bval));
    *buffer=bval;
    break;

  case f_replace:
    bval=*buffer;
    HE(H5LTset_attribute_short(user_attrs.id(),".",name,&bval,nelem));
    break;
  }
}

template <>
void HDF_record_file::declare_attribute(const char* name,short *buffer,
					hsize_t nelem)
{
  switch (file_mode) {
  case f_readonly:
  case f_append:
    HE(H5LTget_attribute_short(user_attrs.id(),".",name,buffer));
    break;

  case f_replace:
    HE(H5LTset_attribute_short(user_attrs.id(),".",name,buffer,nelem));
    break;
  }
}

template <>
void HDF_record_file::declare_attribute(const char* name,int *buffer,
					hsize_t nelem)
{
  switch (file_mode) {
  case f_readonly:
  case f_append:
    HE(H5LTget_attribute_int(user_attrs.id(),".",name,buffer));
    break;

  case f_replace:
    HE(H5LTset_attribute_int(user_attrs.id(),".",name,buffer,nelem));
    break;
  }
}

template <>
void HDF_record_file::declare_attribute(const char* name,long *buffer,
					hsize_t nelem)
{
  switch (file_mode) {
  case f_readonly:
  case f_append:
    HE(H5LTget_attribute_long(user_attrs.id(),".",name,buffer));
    break;

  case f_replace:
    HE(H5LTset_attribute_long(user_attrs.id(),".",name,buffer,nelem));
    break;
  }
}

template <>
void HDF_record_file::declare_attribute(const char* name,long long *buffer,
					hsize_t nelem)
{
  switch (file_mode) {
  case f_readonly:
  case f_append:
    HE(H5LTget_attribute_long_long(user_attrs.id(),".",name,buffer));
    break;

  case f_replace:
    HE(H5LTset_attribute_long_long(user_attrs.id(),".",name,buffer,nelem));
    break;
  }
}
template <>
void HDF_record_file::declare_attribute(const char* name,float *buffer,
					hsize_t nelem)
{
  switch (file_mode) {
  case f_readonly:
  case f_append:
    HE(H5LTget_attribute_float(user_attrs.id(),".",name,buffer));
    break;

  case f_replace:
    HE(H5LTset_attribute_float(user_attrs.id(),".",name,buffer,nelem));
    break;
  }
}
template <>
void HDF_record_file::declare_attribute(const char* name,double *buffer,
					hsize_t nelem)
{
  switch (file_mode) {
  case f_readonly:
  case f_append:
    HE(H5LTget_attribute_double(user_attrs.id(),".",name,buffer));
    break;

  case f_replace:
    HE(H5LTset_attribute_double(user_attrs.id(),".",name,buffer,nelem));
    break;
  }
}

/*****************************************************************************
 *
 * Write
 *
 */

void HDF_record_file::write_header()
{
  for (const field_info &fd : header_fields) {
    HE(H5Dwrite(fd.dataset_id,fd.mem_datatype->id(),H5S_ALL,H5S_ALL,H5P_DEFAULT,
    		fd.buffer));
  }
}

void HDF_record_file::write_record(hsize_t recnum)
{
  // Define the slab (just one point) for i/o
  hsize_t start[1] = {recnum};
  hsize_t count[1] = {1};

  // If recnum>=nrecords, we must resize datasets
  hsize_t newsize[1];
  bool resizing=false;
  if (recnum>=nrecords) {
    nrecords=recnum+1;
    newsize[0]=nrecords;
    resizing=true;
  }
  for (const field_info &fd : record_fields) {
    if (resizing) {
      HE(H5Dset_extent(fd.dataset_id,newsize));
    }
    // First need to select the record through a hyperslab
    // (need a copy of the dataspace)
    hid_t selection;
    HE(selection=H5Dget_space(fd.dataset_id));
    HE(H5Sselect_hyperslab(selection,H5S_SELECT_SET,
    			   start,0,count,0));   
    HE(H5Dwrite(fd.dataset_id,
    		fd.mem_datatype->id(),H5_scalar_s,selection,
    		H5P_DEFAULT,fd.buffer));
    HE(H5Sclose(selection));
  }
}

/*****************************************************************************
 *
 * Read
 *
 */

void HDF_record_file::read_direct(const field_info &fd)
{
  HE(H5Dread(fd.dataset_id,fd.mem_datatype->id(),H5S_ALL,H5S_ALL,H5P_DEFAULT,
	     fd.buffer));
}

void HDF_record_file::read_for_string(const field_info &fd)
{
  char *s;
  HE(H5Dread(fd.dataset_id,fd.mem_datatype->id(),H5S_ALL,H5S_ALL,H5P_DEFAULT,
	     &s));
  *((std::string*)fd.buffer)=s;
  HE(H5free_memory(s));
}

void HDF_record_file::read_header()
{
  for (const field_info &fd : header_fields) {
    // HE(H5Dread(fd.dataset_id,fd.mem_datatype->id(),H5S_ALL,H5S_ALL,H5P_DEFAULT,
    // 		fd.buffer));
    (*(fd.read_function))(fd);
  }
}

void HDF_record_file::read_record(hsize_t recnum)
{
  hsize_t start[1] = {recnum};
  hsize_t count[1] = {1};
  for (const field_info &fd : record_fields) {
    hid_t selection;
    HE(selection=H5Dget_space(fd.dataset_id));
    HE(H5Sselect_hyperslab(selection,H5S_SELECT_SET,
    			   start,0,count,0));   
    HE(H5Dread(fd.dataset_id,fd.mem_datatype->id(),H5_scalar_s,selection,
	       H5P_DEFAULT,fd.buffer));
  }
}


} /* namespace */

