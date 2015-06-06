/*
 * hdf_file.hh -- easier interface for some kinds of HDF5 files
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

#ifndef HDF_FILE_HH
#define HDF_FILE_HH

#include <list>

#include "exception.hh"
#include "cerrors.h"

#include "hdf5.h"
#include "hdf5_hl.h"

#define HE(x)  {if ((x)<0) H5file::handle_HDF_error(HERE);}

namespace glsim {

/******************************************************************************
 *
 * H5file
 *
 */

/** \class H5file
    \ingroup HDF
    \brief Provides error handling and initialization of HDF5 library.

This base class keeps count of the number of open files and calls
`H5open()` before the first one is opened and `H5close()` when the
last one is closed.
*/
class H5file {
public:
  H5file();
  ~H5file();
  hid_t id() const {return hfile_id;}

  virtual void create(const char* fname,unsigned flags,
		      hid_t file_creation_pl=H5P_DEFAULT,
		      hid_t file_access_pl=H5_file_access_p);

  static void handle_HDF_error(Source_context); ///<Prints the HDF error message and thows HDF_error

protected:
  static hid_t H5_scalar_s;  ///< A scalar dataspace
  static hid_t H5_file_access_p;  ///< Suggested file access property list specifying version 1.8

  hid_t hfile_id;

private:
  static int file_count;
} ;

inline void H5file::create(const char* fname,unsigned flags,
			   hid_t file_creation_pl,hid_t file_access_pl)
{
  HE(hfile_id=H5Fcreate(fname,flags,file_creation_pl,file_access_pl));
}


/******************************************************************************
 *
 * HDF_group
 *
 */

/** \class H5group
    \ingroup HDF An HDF5 group

This class holds the handle (`hid`) to a group created or opened with
create() or open() and closes it on destruction.
*/
class H5group {
public:
  H5group() : gid(-1) {}
  ~H5group() {close();}
  hid_t id() const {return gid;}

  void create(hid_t loc_id,const char *name,hid_t link_creation_p=H5P_DEFAULT,
	      hid_t group_creation_p=H5P_DEFAULT,
	      hid_t group_access_p=H5P_DEFAULT)
  {HE(gid=H5Gcreate2(loc_id,name,link_creation_p,group_creation_p,
		     group_access_p));}
  void open(hid_t loc_id,const char *name,hid_t group_access_p=H5P_DEFAULT)
  {HE(gid=H5Gopen2(loc_id,name,group_access_p));}
  void close() {if (gid>=0) HE(H5Gclose(gid)); gid=-1;}
  
private:
  hid_t gid;
} ;

/******************************************************************************
 *
 * H5type
 *
 */

/** \class H5type
    \ingroup HDF

This class is the base of a hierarchy designed to ease the creation of
HDF5 datatypes.  It is intended to hold datatypes defined
independently of a group or file.  The base object holds only the H5
identifier of the type, and closes the type on destruction.

The default constructor creates an invalid type (id=-1) and is only
intended to be called from a descendant class that properly creates
the type.  The public constructor takes a dataset id, which datatype
is copied.

*/
typedef hid_t H5dset;

class H5type {
public:
  H5type(H5dset& ds);  ///<This takes a dataset id and copies the datatype
  ~H5type() {if (tid>=0) HE(H5Tclose(tid));}
  hid_t id() const {return tid;}

protected:
  H5type() : tid(-1) {}

  hid_t tid;

  /// This simply returns the id of the native datatypes as
  /// defined by the HDF5 library.  It is implemented by explicit
  /// instantiation for each type.
  template <typename T> hid_t H5native_datatype();
} ;


inline H5type::H5type(H5dset& ds)
{
  tid=H5Dget_type(ds);
  HE(tid);
}

//
// Explicit instantiations for HDF native datatypes
//

//
// Char and integer
//
template <> inline hid_t H5type::H5native_datatype<char>()
{return H5T_NATIVE_CHAR;}
template <> inline hid_t H5type::H5native_datatype<unsigned char>()
{return H5T_NATIVE_UCHAR;}
template <> inline hid_t H5type::H5native_datatype<short>()
{return H5T_NATIVE_SHORT;}
template <> inline hid_t H5type::H5native_datatype<unsigned short>()
{return H5T_NATIVE_USHORT;}
template <> inline hid_t H5type::H5native_datatype<int>()
{return H5T_NATIVE_INT;}
template <> inline hid_t H5type::H5native_datatype<unsigned int>()
{return H5T_NATIVE_UINT;}
template <> inline hid_t H5type::H5native_datatype<long>()
{return H5T_NATIVE_LONG;}
template <> inline hid_t H5type::H5native_datatype<unsigned long>()
{return H5T_NATIVE_ULONG;}
template <> inline hid_t H5type::H5native_datatype<long long>()
{return H5T_NATIVE_LLONG;}
template <> inline hid_t H5type::H5native_datatype<unsigned long long>()
{return H5T_NATIVE_ULLONG;}

//
// reals
//
template <> inline hid_t H5type::H5native_datatype<float>()
{return H5T_NATIVE_FLOAT;}
template <> inline hid_t H5type::H5native_datatype<double>()
{return H5T_NATIVE_DOUBLE;}
template <> inline hid_t H5type::H5native_datatype<long double>()
{return H5T_NATIVE_LDOUBLE;}

//
// other
//
template <> inline hid_t H5type::H5native_datatype<bool>()
{return H5T_NATIVE_HBOOL;}



/** \class H5simple_type
    \ingroup HDF

This class allows easiy creation of simple types: builtin types, one
and two-dimensional arrays of builtin types, and C++ strings.  Just
provide the type and number of elements upon construction.  For
example:

~~~~~cc
    H5simple_type<double>  T1(N);  // Represents double[N]
    H5simple_type<float,3> T2(N);  // Represents float[N][3]
~~~~~

Note the different syntax for the first dimension, which can be a
variable.  The second, in contrast, must be a constant.  However, you
can make the compiler figure out the constant for you, doing something
like this:

~~~~~cc
     template <typename T,size_t N> inline 
     H5simple_type<T,N> H5simple_type_create(T (&f)[N],hsize_t nel)
     {
       return H5simple_type<T,N>(nel);
     }
~~~~~

*/
template <typename DTYPE,size_t N>
class H5simple_type : public H5type {
public:
  H5simple_type(hsize_t nelements=1);
} ;

template <typename FTYPE,size_t N>
H5simple_type<FTYPE,N>::H5simple_type(hsize_t nelements)
{
  hsize_t dim[2];
  dim[0]=nelements;
  dim[1]=N;
  hid_t native=H5native_datatype<FTYPE>();
  HE(tid=H5Tarray_create2(native,2,dim));
}

template <typename DTYPE>
class H5simple_type<DTYPE,1> : public H5type {
public:
  H5simple_type(hsize_t nelements=1) : H5type()
  {
    hid_t native=H5native_datatype<DTYPE>();
    if (nelements==1) {
      HE(tid=H5Tcopy(native));
    } else {
      hsize_t dim[1];
      dim[0]=nelements;
      HE(tid=H5Tarray_create2(native,1,dim));
    }
  }
} ;

template <>
class H5simple_type<char,1> : public H5type {
public:
  H5simple_type(hsize_t nelements) : H5type()
  {
    HE(tid=H5Tcopy(H5T_C_S1));
    HE(H5Tset_size(tid,nelements));
  }
} ;

template <>
class H5simple_type<std::string,1> : public H5type {
public:
  H5simple_type(hsize_t nelements) {
    HE(tid=H5Tcopy(H5T_C_S1));
    HE(H5Tset_size(tid, H5T_VARIABLE));
  }
} ;

/*****************************************************************************
 *
 * H5space
 *
 */
class H5space {
public:
  H5space() : sid(-1) {}
  ~H5space() {close();}
  hid_t id() {return sid;}

  void close() {if (sid>=0) HE(H5Sclose(sid));}

protected:
  hid_t sid;
} ;

class H5space_simple_base : public H5space {
public:
  H5space_simple_base(int rank) :
    rank_(rank)
  {
    curr_dims=new hsize_t[rank_];
    max_dims=new hsize_t[rank_];
  }
  ~H5space_simple_base() {delete[] curr_dims; delete[] max_dims;}

protected:
  hsize_t rank_,*curr_dims,*max_dims;
} ;

template <int rank>
class H5space_simple : public H5space_simple_base {
public:
  H5space_simple() : H5space_simple_base(rank) {}
  void create(const hsize_t * const cur,const hsize_t * const max);
} ;

template <>
class H5space_simple<3> : public H5space_simple_base {
public:
  H5space_simple() : H5space_simple_base(3) {}
  void create(hsize_t c1,hsize_t c2,hsize_t c3,
	      hsize_t m1,hsize_t m2,hsize_t m3);
} ;

inline void H5space_simple<3>::create(hsize_t c1,hsize_t c2,hsize_t c3,
				      hsize_t m1,hsize_t m2,hsize_t m3)
{
  curr_dims[0]=c1;
  curr_dims[1]=c2;
  curr_dims[2]=c3;
  max_dims[0]=m1;
  max_dims[1]=m2;
  max_dims[2]=m3;
  HE(sid=H5Screate_simple(rank_,curr_dims,max_dims));
}

template <>
class H5space_simple<2> : public H5space_simple_base {
public:
  H5space_simple() : H5space_simple_base(2) {}
  void create(hsize_t c1,hsize_t c2,hsize_t m1,hsize_t m2);
} ;

inline void H5space_simple<2>::create(hsize_t c1,hsize_t c2,
				      hsize_t m1,hsize_t m2)
{
  curr_dims[0]=c1;
  curr_dims[1]=c2;
  max_dims[0]=m1;
  max_dims[1]=m2;
  HE(sid=H5Screate_simple(rank_,curr_dims,max_dims));
}

template <>
class H5space_simple<1> : public H5space_simple_base {
public:
  H5space_simple() : H5space_simple_base(1) {}
  void create(const hsize_t cur,const hsize_t max);
} ;

inline void H5space_simple<1>::create(const hsize_t cur,const hsize_t max)
{
  curr_dims[0]=cur;
  max_dims[0]=max;
  HE(sid=H5Screate_simple(rank_,curr_dims,max_dims));
}

/*****************************************************************************
 *
 * H5set
 *
 */

class H5set {
public:
  H5set() : dset_id(-1) {}
  ~H5set() {close();}
  hid_t id() const {return dset_id;}
  void create(hid_t loc_id,const char *name,hid_t datatype_id,
	      hid_t dataspace_id,hid_t link_creation_pl=H5P_DEFAULT,
	      hid_t dataset_creation_cpl=H5P_DEFAULT,
	      hid_t dataset_access_pl=H5P_DEFAULT);


  void close() {if (dset_id>=0) HE(H5Dclose(dset_id)); dset_id=-1;}
private:
  hid_t  dset_id;
} ;
    
inline void H5set::create(hid_t loc_id,const char *name,hid_t datatype_id,
			  hid_t dataspace_id,hid_t link_creation_pl,
			  hid_t dataset_creation_pl,hid_t dataset_access_pl)
{
  dset_id=H5Dcreate2(loc_id,name,datatype_id,dataspace_id,
		     link_creation_pl,dataset_creation_pl,dataset_access_pl);
  HE(dset_id);
}


/** Exception for HDF5 library errors */
class HDF_error : public glsim::Runtime_error {
public:
  explicit HDF_error(const Source_context& c=Source_context()) :
    Runtime_error("HDF5 error",c) {}
} ;

/******************************************************************************
 *
 * HDF_record_file
 *
 */

/** \class HDF_record_file
    \ingroup HDF

    \brief Main functionality for record-based files under HDF5

This is an abstract class that makes it easy to use simple
header/record files under HDF5.  Simply inherit and declare the fields
belonging to the header and record parts (see declare_header_fields()
and declare_record_fields().  Then the read/write operations can be
used.

Two approaches are possible:

 - flexible buffer organization but unique buffer

 - rigid buffer organization (fixed offsets) but buffer is given at
   read/write time

Right now only the second one is implemented (i.e. field may reside on
arbitrary locations, but on read/write the same buffer is always used.

*/
class HDF_record_file : public H5file {
public:
  /// \name Enums for different options
  /// @{

  /// Open file mode (see open())
  enum mode            {f_readonly, f_replace, f_append};
  /// How to act upon different file versions (see open())
  enum version_require {version_any,version_min,version_exact};
  /// Field can belong to the header or to a record
  enum field_kind      {f_header,f_record};

  /// @}
  /// \name Constructor/destructor and open/close
  /// @{
  HDF_record_file(long version);   ///< Create a file object with given version
  virtual ~HDF_record_file() {close();}
  void               open(const char *fname,const mode& m,
   			  const char *tit="untitled",
			  version_require vr=version_any,long minver=0);
  void               open(hid_t groupid,const mode& m,
   			  const char *tit="untitled",
			  version_require vr=version_any,long minver=0);
  void               close();

  /// @}
  /// \name File information (available while file is open)
  /// @{
  long    file_version() const {return disk_file_version;}
  ///< Version of the opened file
  hsize_t size() const {return nrecords;}
  ///< Size (number of records)
  
  /// @}
  /// \name Read and write
  /// @{
  void write_header();  ///< Write all header fields to disk
  void write_record(hsize_t recnum); ///< Write to specified record number
  void append_record() {write_record(nrecords);} ///< Write new record at the end of the file
  void read_header(); ///< Read all header fields
  void read_record(hsize_t recnum); ///< Read all record fields from record recnum
  /// @}


  /// \class Unsupported_file_version
  /// \ingroup Exceptions
  /// Exception thrown when the file version can't be handled
  class Unsupported_file_version;

protected:

  /// This function must declare al header fields by successive calls
  /// of declare_field().  It is called when opening or creating the
  /// file.  The file opening mode is passed as argument so that
  /// necessary initializations can be done properly when needed.  It
  /// is not mandatory that declaration of header and record fields be
  /// split in two functions; actually both can create any combination
  /// of fields and attributes.  The only difference is that in read
  /// mode, init() calls read_header() after calling
  /// declare_header_fields() so that header information can be used
  /// to declare further fields.  If this information is not needed,
  /// then this function can declare both header and record fields.
  virtual void declare_header_fields(mode fmode)=0;
  /// This must declare al record fields.  It is called when opening
  /// the file.  When this is called all header records are defined,
  /// so that it may call read_header() if header data is needed to
  /// declare the record fields (e.g. to know the size of some array).
  virtual void declare_record_fields(mode fmode)=0;

  /// Call this to declare the fields
  template <typename FTYPE>
  field_kind declare_field(field_kind kind,const char* name,FTYPE *field,
			   hsize_t nelements=1);
  /// This declares an attribute.  Attributes are read or written
  /// immediately after declaration, and as such can serve to store
  /// data needed to setup the rest of the input buffer (e.g. allocate
  /// space for arrays).
  template <typename ATYPE>
  void declare_attribute(const char* name,ATYPE *buffer,hsize_t nelements=1);

private:
  mode         file_mode;
  hsize_t      nrecords;
  hid_t        rootg,H5_frame_s;

  // To create the H5simple_type without explicitly mentioning the type
  // as a template parameter
  template <typename DTYPE> inline H5simple_type<DTYPE,1>*
  H5simple_type_create(DTYPE &f,hsize_t nelements=1)
  {return new H5simple_type<DTYPE,1>(nelements);}

  template <typename DTYPE,size_t N> inline H5simple_type<DTYPE,N>*
  H5simple_type_create(DTYPE (&f)[N],hsize_t nelements=1)
  {return new H5simple_type<DTYPE,N>(nelements);}

private:
  hid_t       file;
  long        file_object_version,disk_file_version;
  std::string title_;
  H5group     user_attrs;

  struct field_info;
  typedef     void (*read_functionT)(const field_info&);
  struct field_info {
    hid_t   dataspace_id;
    H5type  *file_datatype,*mem_datatype;
    hid_t   dataset_id;
    void    *buffer;
    read_functionT read_function;
  } ;
  std::list<field_info> header_fields,record_fields;
  template <typename FTYPE>
  read_functionT read_function(FTYPE *field);

  void        init(const char*,version_require,long);
  static void read_direct(const field_info&);
  static void read_for_string(const field_info&);

} ;

inline HDF_record_file::HDF_record_file(long version) :
  nrecords(0),
  rootg(-1),
  file(-1),
  file_object_version(version)
{
}

class HDF_record_file::Unsupported_file_version : public glsim::Runtime_error {
public:
  explicit Unsupported_file_version(std::string& msg,
				    const Source_context& c=Source_context()) :
    Runtime_error(msg,c) {}
} ;

/**
   This function must be called from declare_header_fields() and
   declare_record_fields() to creates or opens the given field (HDF5
   dataset) depending on the file mode.  Returns the actual kind
   (header or record) of the field opened.

   It supports (through H5simple_type) builtin simple types or one- or
   two- dimensional arrays of builtin types or C++ `std::string`s.  If
   FTYPE is char, then nelements must be the maximum string length.

   \param[in] kind Whether the field is header (f_header) or record
   (f_record).  This only applies when creating the file, when opening
   an existing file the record kind is read from the file.

   \param[in]  name   The (unique) name of the field as will be passed to HDF5

   \param[in] field A pointer to the field (will be used in read and
   write operations)

   \param[in] nelements The number of elements (for a simple array)
   or the maximum string length.  Defaults to 1, meaning a scalar
   (which is different from an array with one element).

   \return The actual kind (f_header or f_record) of the field in the
   opened file (equal to the kind argument if the file is being
   created).

 */
template <typename FTYPE> HDF_record_file::field_kind
HDF_record_file::declare_field(field_kind kind,const char* name,
			       FTYPE *field,hsize_t nelements)
{
  field_info  fd;
  field_kind  actual_kind;
  H5S_class_t dspc;

  switch (file_mode) {

    /* The file exists: must check whether field exists and is
       compatible with request */
  case f_readonly:
  case f_append:
    HE(fd.dataset_id=H5Dopen2(rootg,name,H5P_DEFAULT));
    HE(fd.dataspace_id=H5Dget_space(fd.dataset_id));
    fd.file_datatype=new H5type(fd.dataset_id);
    fd.mem_datatype=H5simple_type_create(*field,nelements);
    dspc=H5Sget_simple_extent_type(fd.dataspace_id);
    switch (dspc) {
    case H5S_SCALAR:
      actual_kind=f_header;
      break;
    case H5S_SIMPLE:
      actual_kind=f_record;
      break;
    default:
      throw glsim::Runtime_error("Unexpected dataspace type",HERE);
      break;
    }
    break;

    /* The file does not exist/is being overwritten: create the field
       from scratch */
  case f_replace:
    fd.file_datatype=H5simple_type_create(*field,nelements);
    fd.mem_datatype=fd.file_datatype;
    switch(kind) {
    case f_header:
      fd.dataspace_id=H5_scalar_s;
      fd.dataset_id=H5Dcreate2(rootg,name,fd.file_datatype->id(),fd.dataspace_id,
			       H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      HE(fd.dataset_id);
      break;

    case f_record:
      fd.dataspace_id=H5_frame_s;
      hsize_t chunks[1]={1};
      hid_t PL=H5Pcreate(H5P_DATASET_CREATE);
      HE(PL);
      HE(H5Pset_chunk(PL,1,chunks));
      fd.dataset_id=H5Dcreate2(rootg,name,fd.file_datatype->id(),fd.dataspace_id,
    			       H5P_DEFAULT,PL,H5P_DEFAULT);
      HE(fd.dataset_id);
      break;
    }
    actual_kind=kind;
    break;
  }

  fd.buffer=field;
  fd.read_function=read_function(field);
  if (actual_kind==f_header) header_fields.push_back(fd);
  else record_fields.push_back(fd);
  return actual_kind;
}

template <typename FTYPE> inline HDF_record_file::read_functionT
HDF_record_file::read_function(FTYPE *field)
{return read_direct;}

template <> inline HDF_record_file::read_functionT
HDF_record_file::read_function<std::string>(std::string *field)
{return read_for_string;}

} /* namespace */

#endif /* HDF_FILE_HH */

