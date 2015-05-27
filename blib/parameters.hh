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

namespace po=boost::program_options;

namespace glsim {

/** \class Parameters
    \ingroup Parameters

To read parameters, one defines a class inherited from Parameters, in
which constructor all needed parameters are declared by calling
parameter_file_options().

To actually read the parameters from the file,
Parameters::parse(char*) must be called, giving it a file name.  The
parsing should be done only once (in the simulation, this is done from
an Environment object).  After parsing, count() and value() may be
called to read parameters as desired.

The backend for parameter reading is Boost::program_options.  Though
we allow that the parameter definition be scattered all over, the
definitios are actually collected in a single static object.  All the
parameters must be defined by the time the parser is called.  As a
result, the library user __must not declare any global Parameters
object,__ or the static member objects may fail to be properly
initialized.

*/
class Parameters {
public:
  static const char* default_scope;

  /// Create the object in the specified scope.
  Parameters(const char* scope_=default_scope) : scope(scope_) {}
  virtual ~Parameters() {}
  /// Parse `parfile` through `Boost::program_options`
  void parse(const char *parfile);
  /// Return the number of values given for the parameter
  int  count(const std::string& s) const;
  /// Get parameter value
  const po::variable_value& value(const std::string& s) const;
  /// Describe all known parameters in the scope
  virtual void show_parameters(std::ostream&) const;
  
protected:
  /// Get an `options_description` object to define parameters
  po::options_description& parameter_file_options();

private:
  typedef std::map<std::string,po::options_description> descmap_t;
  typedef std::map<std::string,po::variables_map> varmap_t;

  std::string      scope;
  static descmap_t description;
  static varmap_t  variables;

  friend class CLParameters;
} ;

/** \class Undefined_parameter
    \ingroup Exception
    \brief Exception thron by Parameters::value() when requested  an undefined parameter
*/
class Undefined_parameter : public glsim::Runtime_error {
public:
  explicit Undefined_parameter(const std::string& param,
			       const Source_context &c=Source_context()) :
  Runtime_error("ERROR: undefined parameter "+param,c) {}
  ~Undefined_parameter() throw() {}
} ;

/** \class Scope_not_parsed
    \ingroup Exception
    \brief Exception thrown by Parameters when a parameter value is requested before calling parse()
*/
class Scope_not_parsed : public glsim::Runtime_error {
public:
  explicit Scope_not_parsed(const std::string& scope,
			       const Source_context &c=Source_context()) :
    Runtime_error("ERROR: requested parameters before scope "+scope+" was parsed",c) {}
  ~Scope_not_parsed() throw() {}
} ;

/** \class Required_parameter_absent
    \ingroup Exception
    \brief Exception thrown requested a required parameter is not defined
*/
class Required_parameter_absent : public glsim::Runtime_error {
public:
  explicit Required_parameter_absent(const std::string& msg,
				     const Source_context &c=Source_context()) :
  Runtime_error("ERROR: "+msg,c) {}
  ~Required_parameter_absent() throw() {}
} ;
/**
To define the parameters, call this function from the children's
constructor to get an `options_description` object from
`Boost::program_options`.  See the example in XX and Boost
documentation for the declaration syntax.
*/
inline po::options_description& Parameters::parameter_file_options()
{
  return description[scope];
}

/**
Returns the number of values given for the specified parameter.  Use
this to determine if an optional parameter has been given (otherwise
value() will throw an exception).
 */
inline int Parameters::count(const std::string& s) const
{
  if (variables.find(scope)==variables.end())
    throw glsim::Scope_not_parsed(scope);
  return variables[scope].count(s);
}


/*****************************************************************************/

/** \class CLParameters
    \ingroup Parameters

CLParametersCL inherits from Parameters, adding functionality to parse
a command line (always relying on boost::program_options).  Only one
CLParameters object should be declared, and it must _not be global._

Command line parameters must be defined by the child's constructor by
calling command_line_options() (used the same way as
parameter_file_options()).  In addition, most programs will define
options called "positional options" in `Boost::program_options`
jargon, i.e. unnamed options that take their value from the position
of the command-line argument instead of using an option syntax (`-o x`
or `--key=value`).  These are defined with positional_options().

This class is abstract because show_usage() is a pure virtual.
Decent command-line parsing requires definition of ``positional
parameters'' (in Boost::program\_options jargon), and a good usage
message to match.  See the test section for an example.
*/
class CLParameters : public Parameters {
public:
  CLParameters(const char* scope=Parameters::default_scope) :
      Parameters(scope)
  {}
  /// Call the Boost parser
  virtual void parse_command_line(int argc,char *argv[],bool merge_options=true);

protected:
  /// Get object to define command line options
  po::options_description&            command_line_options();
  /// Get object to define positional options
  po::positional_options_description& positional_options();
  /// Call Boost::program_options() notify
  void notify();
  /// The (base) name of the program (`argv[0]`), set by parse_commmand_line()
  std::string                               progname;

private:  
  static po::options_description            CLoptions;
  static po::positional_options_description Poptions;
} ;			      

/**
Call this function from the children's constructor to define the
command line parameters through an `options_description` object, as
with parameter_file_options().
*/
inline po::options_description& CLParameters::command_line_options()
{
  return CLoptions;
}

/**
To define positional options call this function and operate through
the returned `po:positional_options_description` object.
*/
inline po::positional_options_description& CLParameters::positional_options()
{
  return Poptions;
}

/**
If parse() was called with merge_options=true, it is necessary to call
Boost::program_options::notify() (which function can throw on errors
like undefined parameters.  If parse() will call notify, but if
parse() is not called after parse_command_line(), then you must call
notify() explicitly.
*/
inline void CLParameters::notify()
{
  po::notify(variables[scope]);
}

/*****************************************************************************/

/** \class SimulationCL
    \ingroup Parameters
    \brief The simulation command line

This class is used to define standard command line for the
simulations distributed with glsim.  It also serves as an additional
example use of CLParameters.  It is meant to use together with an
object of the Enviroment family (i.e. the (unique) StandardCL
object must be created after the Environmnet) and before it is
initialized.  If you want a different syntax for the command line you
can define another class to replace this, just be sure to add the
parameters defined in the constructor, either as command-line or
parameter-file options.
*/
class SimulationCL : public CLParameters {
public:
  SimulationCL(const char* name,const char* copyright,const char *scope=Parameters::default_scope);
  void parse_command_line(int argc,char *argv[]);
  virtual void show_usage();

private:
  std::string name_and_ver,copyr;

  void show_version();
} ;

/// In addition to Early_stop_required, ParametersCL can throw Usage_error.
class Usage_error : public glsim::Runtime_error {
public:
  explicit Usage_error() :
  Runtime_error("usage error, try -h or --help",Source_context())
  {}
} ;

/*****************************************************************************/

/** \class UtilityCL
    \ingroup Parameters
    \brief Common command-line options for utility programs

This provides a parser and defines options common to all utility or
analysis programs distributed with glsim.
*/
class UtilityCL : public CLParameters {
public:
  UtilityCL(const char* name,const char *scope=Parameters::default_scope);
  virtual void parse_command_line(int argc,char *argv[]);

private:
  std::string name_and_ver;

  void show_version();
  virtual void show_usage()=0;
} ;

/** \brief Catch exceptions expected by standard utilities

This simply calls the given function wmain, wrapping the call in a
try...catch block that catches the glsim library exceptions and prints
the backtrace.
 */
int UtilityEC(int argc,char *argv[],void (*wmain)(int,char**));

} /* namespace */

#endif /* PARAMETERS_HH */
