/*
 * scontext.cc -- source context methods
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

#include "scontext.hh"

/*

* \file scontext.cc

    Provides overload of operator<<

 */


/** \namespace glsim  Namespace for all glsim classes and functions
 */
namespace glsim {

/** \ingroup Error

 Overload of `operator<<` to output bactrace to a stream

*/
std::ostream& operator<<(std::ostream& o,const glsim::Backtrace& b)
{
  char **strings = backtrace_symbols(b.cbuffer, b.nptrs);
  o << "===== BACKTRACE START ==========\n";
  if (strings == 0)
    o << "ERROR: Failed to translate backtrace\n";
  else
    for (int j = 0; j < b.nptrs; j++)
      o << strings[j] << '\n';
  o << "===== BACKTRACE END ==========\n";
  free(strings);
  return o;
}

}

