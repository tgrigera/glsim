/*
 * parameters.hh -- declarations for parameters classes
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

#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#include <string>
#include <iostream>
#include <map>
#include <boost/program_options.hpp>

#include "exception.hh"
#include "cerrors.h"

namespace glsim {

namespace po=boost::program_options;

/** \class Parameters
    \ingroup Simulation

To read parameters, the user declares a class inherited from
[[Parameters]].  The constructor of the derived class must declare
the parameters to be read by calling [[parm_file_options]], which is
an object of type [[options_description]] from
[[boost::program_options]].  See the example in [[test]] and the Boost
documentation for the declaration syntax.  Optionally, a scope can be
given when the object is created.  It is legal to define two objects
of the same class with different scopes.

To actually read the parameters from the file, one must call
[[Parameters::parse(char*)]] passing it a file name.  The parsing
should be done only once (in the simulation, from an [[Environent]]
object, which see).  After that, [[Parameters::count]] and
[[Parameters::value]] may be called to load parameters as desired.

The backend for parameter reading is Boost::program\_options.  Though
we allow that the parameter definition be scattered all over, the
definitios are actually collected in a single static object.  All the
parameters must be defined by the time the parser is called.  As a
result, the library user \emph{must not declary any global
  [[Parameter]] object,} or the static member objects may fail to be
properly initialized.

*/
class Parameters {
public:
  static const char* default_scope;

  Parameters(const char* scope_=default_scope) : scope(scope_) {}
  void parse(const char *parfile);
  int  count(const std::string& s) const;
  const po::variable_value& value(const std::string& s) const;
  virtual void show_parameters(std::ostream&) const;
  
protected:
  po::options_description& parm_file_options();

private:
  typedef std::map<std::string,po::options_description> descmap_t;
  typedef std::map<std::string,po::variables_map> varmap_t;

  std::string      scope;
  static descmap_t description;
  static varmap_t  variables;

  friend class ParametersCL;
} ;

/// [[Parameters::value()]] will throw an/ exception when an undefined parameter is requested.
class Undefined_parameter : public glsim::Runtime_error {
public:
  explicit Undefined_parameter(const std::string& param,
			       const Source_context &c=Source_context()) :
  Runtime_error("ERROR: undefined parameter "+param,c) {}
  ~Undefined_parameter() throw() {}
} ;

class Scope_not_parsed : public glsim::Runtime_error {
public:
  explicit Scope_not_parsed(const std::string& scope,
			       const Source_context &c=Source_context()) :
    Runtime_error("ERROR: requested parameters before scope "+scope+" was parsed",c) {}
  ~Scope_not_parsed() throw() {}
} ;

/// To define the options, call this from the child's constructor to get
/// an [[options_description]] object from Boost.
inline po::options_description& Parameters::parm_file_options()
{
  return description[scope];
}

inline int Parameters::count(const std::string& s) const
{
  if (variables.find(scope)==variables.end())
    throw glsim::Scope_not_parsed(scope);
  return variables[scope].count(s);
}


/** \class ParametersCL
    \ingroup Simulation
[[ParametersCL]] inherits from [[Parameters]], adding functionality to
parse a command line, always relying on [[boost::program_options]]'
facilities.  Only one [[ParametersCL]] object should be declared, and
it must \emph{not} be global. %'

This class is abstract because [[show_usage]] is a pure virtual.
Decent command-line parsing requires definition of ``positional
parameters'' (in Boost::program\_options jargon), and a good usage
message to match.  See the test section for an example.
*/
class ParametersCL : public Parameters {
public:
  ParametersCL(const char* scope=Parameters::default_scope);
  virtual void parse_command_line(int argc,char *argv[],bool require_parameter_file=true,
				  bool use_parameter_file=true);
  virtual void show_usage()=0;

protected:
  static po::options_description            command_line_options;
  static po::positional_options_description pos;
  std::string                               progname;
} ;			      


/** \class StandardCL
    \ingroup Simulation
    \brief The standard command line

The following class is used to define standard command line for the
simulations distributed with glsim.  It also serves as an additional
example use of ParametersCL.  It is meant to use together with an
object of the Enviroment family (i.e. the (unique) StandardCL
object must be created after the Environmnet) and before it is
initialized.  If you want a different syntax for the command line you
can define another class to replace this, just be sure to add the
parameters defined in the constructor, either as command-line or
parameter-file options.  The parameters corresponding to the options
`-c`, `-i`, and `-f` below are only used in the STANDARD INIT???.
*/
class StandardCL : public ParametersCL {
public:
  StandardCL(const char *scope=Parameters::default_scope);
  void show_usage();
} ;

/// In addition to Early_stop_required, ParametersCL can throw Usage_error.
class Usage_error : public glsim::Runtime_error {
public:
  explicit Usage_error() :
  Runtime_error("usage error, try -h or --help",Source_context())
  {}
} ;


} /* namespace */

#endif /* PARAMETERS_HH */
