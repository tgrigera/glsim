/*
 * scontext.hh -- source context
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

#ifndef SCONTEXT_HH
#define SCONTEXT_HH

#include <execinfo.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>

namespace glsim {



/** \class Backtrace
    \ingroup Scontext
    \brief Store a stack trace

This class stores a backtrace, which can later be printed to a stream
using the overloaded insertion operator.  Note that compiling with
debugging options and requesting more detailed symbol tables in the
executable (e.g. with `-g` and `-rsymbol` in `gcc`) and will usually
give richer backtrace information (but in any case try processing with
addr2line with the aid of the script below).

## Interpreting the backtrace.

The script `readbt.sh` extracts the backtrace produced by operator<< 
and processes it with `addr2line` to translate the return addresses to
source lines and function names.  It can be run like this:

     ./progname 2>&1 | readbt.sh ./progname

(note that you must give `progname` again as an argument to
`readbt.sh`).  Alternatively you can copy the programs output to a
file and redirect input to `readbt.sh` from that file.

*/
class Backtrace {
public:
  Backtrace();  ///<Stores a backtrace on construction

  friend std::ostream& operator<<(std::ostream&,const Backtrace&);

private:
  int              nptrs;
  static const int cbuffer_len = 100;
  void             *cbuffer[cbuffer_len];
} ;

inline Backtrace::Backtrace()
{
  nptrs = backtrace(cbuffer,cbuffer_len);
}
  
/// \ingroup Scontext
#define HERE glsim::Source_context(__FILE__,__LINE__,__func__)

/** \class Source_context
    \ingroup Scontext
    \brief Store source filename, line number and function name

This class holds source filename, line and function information stored
as a string.  It is intended mainly to be used with the `HERE` macro
to store local context (to report an error or throw an exception).  It
stores automatically a backtrace than can be requested through the
backtrace() method.  This feature is used by the exceptions defined
below.

*/
class Source_context {
private:
  std::string desc_;
  Backtrace    bt;
  
public:
  Source_context(const char *FILE=0,int LINE=0,const char *func=0);
  const std::string& description() const {return desc_;}
  const Backtrace& backtrace() const {return bt;}
} ;
  
inline Source_context::Source_context(const char *FILE,int LINE,const char *func)
{
  if (FILE==0) {
    desc_="";
    return;
  }

  std::stringstream ss;
  ss << FILE << ':' << LINE;
  if (func) ss << " (" << func << ')';
  ss << ": ";
  desc_ = ss.str();
}

inline std::ostream& operator<<(std::ostream& o,const Source_context& c)
{
  o << c.description();
  return o;
}

}

#endif /* SCONTEXT_HH */
