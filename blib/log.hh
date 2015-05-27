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

/** The four verbosity levels */
enum loglevel {debug=0, info=1, warn=2, error=3};

/** \class logger
    \ingroup Logging

Only one global object of this class should exist, defined in the
library and declared in this header.

The output is done through operator():

     glsim::logs(glsim::warn) << "WARNING\n";

*/
class logger {
public:
  logger();  ///< Creates the log stream defaulting to std::cout at debug level
  std::ostream& operator()(loglevel) const;  ///<Return reference to the output stream at the given level

  logger& set_stream(std::ostream&,loglevel);  ///< Chose the log stream and verbosity
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


/**
The constructor initializes the logstream at the highest level and
outputing to `std::cout`, and sets up the null stream.  The easiest
way to do the last appears to be the trick suggested in StackOverflow
(http://stackoverflow.com/questions/6240950/platform-independent-dev-null-in-c/6240980\#6240980):
define a `std::ostream` initialized with a null pointer.  This creates
a stream without a `streambuf` buffer, so that the stream is left in
an error state and never outputs anything (I also believe that it
skips all formatting, thus calls to the insertion operator should be
rather fast).
*/
inline logger::logger() :
  nullstream(0)
{
  set_stream(std::cout,debug);
}

/** 
operator() simply returns a reference to the actual stream that will
do the output.  This stream can be

 1. an actual output stream such as `std::cerr` or a `std::ofstream` linked to a file,

 2. a tee stream to implement writing to two sinks with just one
    insertion operator, or

 3. a null stream that discards everything sent to it.
*/
inline std::ostream& logger::operator()(loglevel l) const
{
  return *(logstream[l]);
}
  
} /* namespace */

#endif /* LOG_HH */
