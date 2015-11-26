/*
 * exception.hh -- exception hierarchy for glsim
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

#ifndef EXCEPTION_HH
#define EXCEPTION_HH

#include <sstream>
#include <stdexcept>

#include <errno.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include "scontext.hh"

namespace glsim {

  
/** \class Logic_error
    \ingroup Exceptions

    \brief A std::logic_error descendant that stores source context
*/
class Logic_error : public std::logic_error {
private:
  const Source_context scontext;

public:
  Logic_error(const std::string& desc,const Source_context &c=Source_context() ) :
    std::logic_error(c.description()+desc),
    scontext(c)
  {}
  const Backtrace &backtrace() const {return scontext.backtrace();}
  ~Logic_error() throw () {}
} ;

/** \class Runtime_error
    \ingroup Exceptions

    \brief A std::runtime_error descendant that stores source context
*/
class Runtime_error : public std::runtime_error {
private:
  const Source_context scontext;

public:
  Runtime_error(const std::string& desc,const Source_context &c=Source_context() ) :
     std::runtime_error(c.description()+desc),
     scontext(c)
  {}
  const Backtrace &backtrace() const {return scontext.backtrace();}
  ~Runtime_error() throw () {}
} ;


/** \group logic_error_children Exceptions deriving from Logic_error
    \ingroup Exceptions

\class Unimplemented
\brief Signal unimplemented feature.
\ingroup logic_error
*/
class Unimplemented : public Logic_error {
public:
  explicit Unimplemented(const std::string& feature="unkown",
			 const Source_context &c=Source_context() ) :
    Logic_error("unimplemented: "+feature,c)
    {}
} ;

/** \class Invalid_value
    \ingroup Exceptions
 */
class Invalid_value : public Logic_error {
public:
  explicit Invalid_value(const std::string& value,const std::string& par,
			 const Source_context &c=Source_context() ) :
    Logic_error("invalid value ("+value+") for parameter "+par,c)
  {}
  ~Invalid_value() throw () {}
} ;


/** \class Out_of_range
    \ingroup Exceptions
*/
class Out_of_range : public Logic_error {
public:
  explicit Out_of_range(const Source_context &c=Source_context() ) :
    Logic_error("out of range",c) {}
  ~Out_of_range() throw () {}
} ;


/** Runtime_error children
 \class Invalid_operation
 \ingroup Exceptions
*/
class Invalid_operation : public Runtime_error {
public:
  explicit Invalid_operation(const std::string& operation=std::string(),
			     const Source_context &c=Source_context()) :
    Runtime_error("invalid operation" + (operation.empty() ? "" : ": " + operation),
		  c)
  {}
} ;

/** Runtime_error children
 \class Open_file_error
 \ingroup Exceptions
*/
class Open_file_error : public Runtime_error {
public:
  explicit Open_file_error(const std::string& file,
			   const Source_context &c=Source_context()) :
    Runtime_error("Cannot open "+file+": "+strerror(errno) )
  {}
} ;

/** Runtime_error children
    \class Clib_error
    \ingroup Exceptions
*/
class Clib_error : public Runtime_error {
public:
  explicit Clib_error(const Source_context &c=Source_context()) :
    Runtime_error(std::string("C library error: ") + strerror(errno), c)
    {}
} ;

/** Runtime_error children
    \class GSL_error
    \ingroup Exceptions
*/
class GSL_error : public Runtime_error {
private:
  std::string m;
public:
  explicit GSL_error(const int gsl_ecode,
		     const Source_context &c=Source_context()) : 
    Runtime_error(std::string("GSL error: ") + gsl_strerror(gsl_ecode))
  {}
  ~GSL_error() throw() {}
} ;


/** \class Early_stop
    \ingroup Exceptions

This is not actually an error, but a condition that requires early
stop.  It should be caught but not reported (i.e.\ what() not called).

*/
class Early_stop : public Runtime_error {
public:
  explicit Early_stop() :
    Runtime_error("Early stop",Source_context())
  {}
} ;



} /* namespace */


#endif /* EXCEPTION_HH */
