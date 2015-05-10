/*
 * log.hh -- a simple logging stream
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

#ifndef LOG_HH
#define LOG_HH

#include <iostream>

namespace glsim {

enum loglevel {debug=0, info=1, warn=2, error=3};

/** \class logger
    \ingroup Logging
*/
class logger {
public:
  logger();
  std::ostream& operator()(loglevel) const;

  logger& set_stream(std::ostream&,loglevel);
  logger& set_additional_stream(std::ostream&,loglevel);

private:
  struct {std::ostream* str; loglevel level;} stream1,stream2;
  std::ostream* logstream[4];
  std::ostream nullstream;
} ;

#ifdef _LOGGING_CC_
logger logs;
#else
extern logger logs;
#endif

inline logger::logger() :
  nullstream(0)
{
  set_stream(std::cout,debug);
}

inline std::ostream& logger::operator()(loglevel l) const
{
  return *(logstream[l]);
}
  
} /* namespace */

#endif /* LOG_HH */
