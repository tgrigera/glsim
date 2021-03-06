/*
 * hdf_file.cc -- class interface for HDF5 files (definitions)
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
    H5Pset_libver_bounds(H5_file_access_p, H5F_LIBVER_EARLIEST, H5F_LIBVER_LATEST);
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

  // Create the frame (unlimited) dataspace
  hsize_t cur[1]={0};
  hsize_t max[1]={H5S_UNLIMITED};
  H5_frame_s=H5Screate_simple(1,cur,max);
  HE(H5_frame_s);

  // Open and check version
  switch(file_mode) {
  case f_readonly:
  case f_append:
    HE(H5_maj=H5Aopen(rootg,"structure_version",H5P_DEFAULT));
    HE(H5Aread(H5_maj,H5T_NATIVE_LONG,&disk_file_version));
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

  //  HE(H5Dvlen_reclaim(fd.mem_datatype->id(),space,H5P_DEFAULT,s));   // Some form of this functionmay be needed to reclaim space when using variable length strings, but documentation is obscure and I don't know what dataspace has to be specified (H5S_ALL is no good)

  HE(H5free_memory(s));
}

void HDF_record_file::read_header()
{
  for (const field_info &fd : header_fields) {
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
    HE(H5Sclose(selection));
  }
}


/******************************************************************************
 *
 * H5_multi_file
 *
 */

H5_multi_file::H5_multi_file(std::vector<std::string> filelist,HDF_record_file &fileob) :
  own_ptr(false),
  filep(&fileob)
{
  currpos.file=0;
  currpos.file_rec=0;
  currpos.global_rec=0;
  filedesc.resize(filelist.size());
  totrecs=0;
  for (int i=0; i<filelist.size(); ++i) {
    filedesc[i].name=filelist[i];
    filep->open_ro(filedesc[i].name.c_str());
    filedesc[i].Nrec=filep->size();
    if (filedesc[i].Nrec==0) filedesc[i].Nrec=1;
    filedesc[i].first_rec=totrecs;
    totrecs+=filedesc[i].Nrec;
    filep->close();
  }
  open(0);
}

void H5_multi_file::open(int filen)
{
  currpos.file=filen;
  currpos.file_rec=0;
  currpos.global_rec=filedesc[filen].first_rec;
  filep->close();
  filep->open_ro(filedesc[currpos.file].name.c_str());
  filep->read_header();
}

bool H5_multi_file::open_next()
{
  if (currpos.file==filedesc.size()-1) return false;
  open(currpos.file+1);
  return true;
}

// note that there may be files with pure header and no record
// data, but these must be counted has having one record anyway
// because the header data may change from file to file, and
// fields that are record fields in some files may be header fields
// in other files
bool H5_multi_file::read()
{
  if (eof()) return false;
  if (currpos.file_rec<filep->size())
    filep->read_record(currpos.file_rec);
  ++currpos.file_rec;
  ++currpos.global_rec;
  return true;
}

bool H5_multi_file::eof()
{
  if (currpos.file_rec<filep->size() || currpos.file_rec==0) return false;
  return !open_next();
}

H5_multi_file& H5_multi_file::seek(hsize_t rec)
{
  if (rec>totrecs) rec=totrecs;
  if (rec<0) rec=0;

  currpos.global_rec=rec;
  int newfile=filedesc.size()-1;
  while (rec<filedesc[newfile].first_rec) --newfile;
  currpos.file_rec=rec-filedesc[newfile].first_rec;
  if (currpos.file!=newfile) {
    currpos.file=newfile;
    filep->close();
    filep->open_ro(filedesc[currpos.file].name.c_str());
    filep->read_header();
  }
  return *this;
}

H5_multi_file& H5_multi_file::rewind()
{
  if (currpos.file==0) {
    currpos.file_rec=0;
    currpos.global_rec=0;
  } else {
    filep->close();
    open(0);
  }
  return *this;
}

} /* namespace */

