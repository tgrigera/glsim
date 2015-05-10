/*
 * parameters.cc -- definitions for parameters classes
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


#include <fstream>

#include "log.hh"
#include "parameters.hh"

namespace glsim {
  
// This is where the descriptions and values are stored.  To implement
// scopes we simply use a [[std::map]].

const char* Parameters::default_scope="[default]";
Parameters::descmap_t Parameters::description;
Parameters::varmap_t  Parameters::variables;

po::options_description            ParametersCL::command_line_options;
po::positional_options_description ParametersCL::pos;

// The parse function just properly calls the Boost component.

void Parameters::parse(const char *parfile)
{
  po::store(po::parse_config_file<char>(parfile,parm_file_options()),
	    variables[scope]);
  po::notify(variables[scope]);
  glsim::logs(glsim::info) << "Parameters read from " << parfile << "\n";
}
  
// With the next methods one can retrieve the values of the
// options. See the Boost documentation for usage of the
// [[po::variable_value]] object.

const po::variable_value& Parameters::value(const std::string& s) const 
{
  if (variables.find(scope)==variables.end())
    throw Scope_not_parsed(scope);
  if (variables[scope].count(s)>0)
    return variables[scope][s];
  else
    throw Undefined_parameter(s);
}

// A help screen can be printed easily thanks to
// [[boost::program_options]]' facilities.  Override if you prefer a
// custom message.

void Parameters::show_parameters(std::ostream& o) const
{
  o << description[scope];
}




// The constructor declares three command-line options ([[--help]] or
// [[-h]] to request usage help, [[--parameter-help]] to show all the
// accepted parameters, and the optional [[parameter_file]] to name the
// [[.ini]] file to be used).  The rest must be added by the user through
// the [[command_line_options]] object.  Positional options can be
// declared through the [[pos]] protected object, as documented in the
// Boost library (see also the examples in the test section).

ParametersCL::ParametersCL(const char* scope) :
  Parameters(scope)
{
  command_line_options.add_options()
    ("help,h",po::bool_switch(),"help with usage")
    ("parameter-help",po::bool_switch(),"show accepted parameters")
    ("parameter_file",po::value<std::string>(),"specifiy parameter (.ini) file")
    ;
}

// This method parses the command line and acts upon the help options,
// displaying the requested help and throwing an exception.  Also, if
// [[use_parameter_file]] is true and the [[control_file]] option is
// given, the named file is parsed.  With [[require_parameter_file]], an
// exception will be thrown if the parameter file is not given.  If this
// behavior is not wanted, the method can be overriden.  To check for the
// legality of the command line, in simple cases it will suffice to mark
// some parameters as [[required()]].  In more complex cases a
// [[parse_command_line]] can be written that calls this one to parse the
// command line and then checks that all required command-line parameters
// have been read and are consistent.

// This method stores the program name in [[progname]] for the benefit of
// [[show_usage]], if you override it be sure to initialize it also.

void ParametersCL::parse_command_line(int argc,char *argv[],bool require_parameter_file,
				      bool use_parameter_file)
{
  progname=basename(argv[0]);
  std::string parameter_file;

  try {

    command_line_options.add(parm_file_options());

    po::store(po::command_line_parser(argc,argv).options(command_line_options).
              positional(pos).run(),variables[scope]);
  
    if (value("help").as<bool>()) {
      show_usage();
      throw Early_stop();
    }
    if (value("parameter-help").as<bool>()) {
      show_parameters(std::cerr);
      throw Early_stop();
    }
    if (count("parameter_file")>0) {
      if (require_parameter_file || use_parameter_file)
	parameter_file=value("parameter_file").as<std::string>();
	parse(parameter_file.c_str());
    } else {
      if (!require_parameter_file) po::notify(variables[scope]);
      else throw Usage_error();
    }

  } catch (po::too_many_positional_options_error& e) {
    throw Usage_error();
  } catch (po::invalid_command_line_syntax& e) {
    throw Usage_error();
  } catch (po::invalid_command_line_style& e) {
    throw Usage_error();
  } catch (po::required_option& e) {
    throw Usage_error();
  } catch (po::reading_file& e) {
    throw Runtime_error("Cannot read parameter file "+parameter_file);
  }
}



//<<Standard command line methods>>=
StandardCL::StandardCL(const char* scope) :
  ParametersCL(scope)
{
  command_line_options.add_options()
    ("initial_infix",po::value<std::string>()->required())
    ("final_infix",po::value<std::string>()->required())
    ("configuration-init,c",po::value<std::string>())
    ("ignore-partial-run,i",po::bool_switch())
    ("force-overwrite,f",po::bool_switch())
    ;
  pos.add("parameter_file",1).add("initial_infix",1).add("final_infix",1);
}

void StandardCL::show_usage()
{
  std::cerr << "usage: " << progname << " [options] parameter_file initial_infix final_infix\n\n"
    << "parameter_file is an ASCII file (.ini style) with the definition of the\n"
    << "simulation parameters.  initial_infix and final_infix will be used to form\n"
    << "the names of the input and output files (check the help for parameters\n"
    << "configuration_file_prefix) unless overridden by setting\n"
    << "initial_configuration_file, etc explicitly.\n\n"
    << "When starting from scratch, give +++ as initial_inifx to indicate that\n"
    << "configuration and environment must not be read from an earlier run.\n"
    << "In this case you may wish to give option -c below.\n\n"
    << "Options:\n"
    << "  -c,--configuration-init string\n"
    << "                          Will pass \"string\" to the Configuration\n"
    << "                          init method, to request inizialization of a default\n"
    << "                          configuration.  The string may be interpreted as\n"
    << "                          a filename.\n"
    << "  -f,--force-overwrite    Overwrite any eventual files with completed runs\n"
    << "  -i,--ignore-partial-run Ignore eventual files with partially completed runs\n"
    << "  -h,--help               Show this help\n"
    << "  --parameter-help        Show accepted parameters\n"
	 ;
}






} /* namespace */

