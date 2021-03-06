/*
 * mfile.hh -- Treat multiple files as one (for input)
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use glsim to produced published work, or if you redistribute a
 * modified version of glsim, or code based on glsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * glsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding.  * glsim
 * distribution.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */


#ifndef MFILE_HH
#define MFILE_HH

#include <string>
#include <vector>
#include <cstdio>

#include "exception.hh"

namespace glsim {

/** \class MFILE
    \ingroup FILE
    \brief Manipulate multiple files as one large file (for input)

    This class gives random, read-only, access to multiple files,
    presenting them as one large file.  This uses C-style i/o.

    Reading is done with the C-library functions, since MFILE can be
    converted to a FILE*, but files are never opened or closed
    explicitly.  Before reading, check for end of file with eof():
    this returns true at the end of the last file.  Otherwise it
    returns false, but if it is at the end of some other file, it
    automatically opens the next one.  It doesn't need to read past
    eof.

    For example, a simple 'cat' would be something like

    ~~~~~{.cc}
    MFILE mf(flist);
    while (!mf.eof()) {
      char buff[200];
      fgets(buff,200,mf);
      std::cout << buff;
    }
    ~~~~~

    You don't need a second read at the end to trigger eof().

 */
class MFILE {
public:
  MFILE(std::vector<std::string> filelist);  ///< Creates MFILE and opens first file.  If filelist is empty, stdin is used
  ~MFILE() noexcept(false);  ///< Close all files and destroy

  bool eof();   ///< Returns true when last character of last file has been read.  Does not need past-eof-read to detect EOF
  void rewind(); ///< Return to the beginning of the first file
  void set_mark(int m=0); ///< Save the current position
  void goto_mark(int m=0); ///< Go to the saved mark
  operator FILE*() ///< Cast so that MFILE can be used in place of a FILE*
  {return filep;}

private:
  std::vector<std::string> filelist;
  unsigned int fileno;
  FILE     *filep;
  int      mark_files;
  long     mark_fpos;

  void open_first();
  void open_next();
} ;

  inline MFILE::~MFILE() noexcept(false)
{
  if (filep)
    if (fclose(filep)) throw Clib_error(HERE);
}

} /* namespace */

#endif /* MFILE_HH */

