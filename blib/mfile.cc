/*
 * mfile.cc -- Treat multiple files as one (for input)
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

#include "mfile.hh"

namespace glsim {

MFILE::MFILE(std::vector<std::string> filelist_) :
  filelist(filelist_),
  fileno(-1),
  filep(0)
{
  if (filelist.empty()) filep=stdin;
  else open_first();
}

void MFILE::open_first()
{
  fileno=0;
  if ( (filep=fopen(filelist[fileno].c_str(),"r"))==0) throw Clib_error(HERE);
}

void MFILE::open_next()
{
  if (filep==stdin) return;
  if (fileno+1<filelist.size()) {
    if (fclose(filep)) throw Clib_error(HERE);
    fileno++;
    if ( (filep=fopen(filelist[fileno].c_str(),"r"))==0) throw Clib_error(HERE);
  }
}

void MFILE::rewind()
{
  if (fileno==0) {::rewind(filep); return;}
  if (fclose(filep)) throw Clib_error(HERE);
  open_first();
}

bool MFILE::eof()
{
  if (ungetc(fgetc(filep),filep)!=EOF) return false;
  open_next();
  return ungetc(fgetc(filep),filep)==EOF;
}

void MFILE::set_mark(int m)
{
  if (m!=0) throw Unimplemented("Multiple marks in MFILE");
  mark_files=fileno;
  mark_fpos=ftell(filep);
}

void MFILE::goto_mark(int m)
{
  if (m!=0) throw Unimplemented("Multiple marks in MFILE");
  if (fileno!=mark_files) {
    fileno=mark_files-1;
    open_next();
  }
  fseek(filep,mark_fpos,SEEK_SET);
  eof();
}

}
