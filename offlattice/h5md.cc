/*
 * h5md.cc -- class interface for H5MD trajectory files (definitions)
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

#include "../config.h"

#include "h5md.hh"

namespace glsim {

H5MD::H5MD(const char *fname,const char *author_name,const char *author_email)
{
  if (author_name!=0) {  // create mode
    HE(hfile_id=H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5_file_access_p));
    h5md.create(hfile_id,"h5md");
    int v[2]={1,0};
    HE(H5LTset_attribute_int(h5md.id(),".","version",v,2));
    author.create(h5md.id(),"author");
    HE(H5LTset_attribute_string(author.id(),".","name",author_name));
    if (author_email!=0)
      HE(H5LTset_attribute_string(author.id(),".","email",author_email));
    creator.create(h5md.id(),"creator");
    HE(H5LTset_attribute_string(creator.id(),".","name","glsim"));
    HE(H5LTset_attribute_string(creator.id(),".","version",PACKAGE_VERSION));
  } else {
    HE(hfile_id=H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT));
    h5md.open(hfile_id,"h5md");
  }

}

} /* namespace */
