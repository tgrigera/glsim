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
#include "config.h"

namespace glsim {
  
/*****************************************************************************
 *
 * class Parameters
 *
 */

// This is where the descriptions and values are stored.  To implement
// scopes we simply use a std::map.

const char*           Parameters::default_scope="[default]";
Parameters::descmap_t Parameters::description;
Parameters::varmap_t  Parameters::variables;

// The parse function just properly calls the Boost component.

void Parameters::parse(const char *parfile)
{
  po::store(po::parse_config_file<char>(parfile,parameter_file_options()),
	    variables[scope]);
  po::notify(variables[scope]);
  glsim::logs(glsim::info) << "Parameters read from " << parfile << "\n";
}
  
/**
  This returns a reference to an appropriate
  `Boost::program_options::variable_value` object, which can be used
  to retrieve the value of the parameter (see the
  `Boost::program_options` documentation for usage of the
  `variable_value` object.
*/
const po::variable_value& Parameters::value(const std::string& s) const 
{
  if (variables.find(scope)==variables.end())
    throw Scope_not_parsed(scope);
  if (variables[scope].count(s)>0)
    return variables[scope][s];
  else
    throw Undefined_parameter(s);
}

/**
   This prints a description of all the parameters defined in the
   scope through `Boost::program_options`.  Override if you prefer a
   custom message.
*/
void Parameters::show_parameters(std::ostream& o) const
{
  o << description[scope];
}

/*****************************************************************************
 *
 * class CLParameters
 *
 */

po::options_description            CLParameters::CLoptions;
po::positional_options_description CLParameters::Poptions;
/**

Parse the command line and set the program name in progname.
Optionally add all parameter file options to the command line options:
in this way parameters can be given in the command line.  If the
command line is parsed before any parameter files, the command-line
parameters will override those given in the parameter file.

\param[in] argc,argv     Argument count and values as passed to main()

\param[in] merge_options If `true`, add the parameter file options to
                         the command line options.  In this case, the
                         Boost notify function is not called, since it
                         may still be necessary to parse a parameter
                         file (notfiy must be called only when all
                         sources have been parsed, since it can throw
                         exceptions related to missing options).
                         Notify will then be called by parse().

*/
void CLParameters::parse_command_line(int argc,char *argv[],
				      bool merge_options)
{
  progname=basename(argv[0]);
  if (merge_options)
    CLoptions.add(parameter_file_options());

  po::store(po::command_line_parser(argc,argv).options(CLoptions).
	    positional(Poptions).run(),variables[scope]);
  if (!merge_options) po::notify(variables[scope]);
}

/*****************************************************************************
 *
 * class SimulationCL
 *
 */

/**
The constructor declares six command-line options (`--help` or `-h`,
`--version`, `--list-parameters`, `-c`, `-i`, and `-f`) plus three
positional options for the parameter file and initial and final infix.
More can be added by the user through command_line_options().

 \param[in] name_and_version  Name and version string to print with --version
 \param[in] copyright         Copyright string to print with --version
*/
SimulationCL::SimulationCL(const char* name_and_version,const char* copyright,const char* scope) :
  CLParameters(scope),
  name_and_ver(name_and_version),
  copyr(copyright)
{
  command_line_options().add_options()
    ("help,h",po::bool_switch(),"help with usage")
    ("version",po::bool_switch(),"print version and exit")
    ("list-parameters",po::bool_switch(),"show accepted parameters")
    ("parameter_file",po::value<std::string>(),"read further options/parameters from (.ini) file")
    ("initial_infix",po::value<std::string>()->required())
    ("final_infix",po::value<std::string>()->required())
    ("configuration-init,c",po::value<std::string>())
    ("ignore-partial-run,i",po::bool_switch())
    ("force-overwrite,f",po::bool_switch())
    ;
  positional_options().add("parameter_file",1).add("initial_infix",1).
    add("final_infix",1);
}

/**
This method parses the command line (calling
CLParameters::parse_comand_line()) and acts upon the help options
(`--help`, `--version`, and `--list-parameters`) displaying the
requested help and throwing an exception to request an early stop.
Also, if the `parameter-file` option is given, the named file is
parsed.  When `require_parameter_file` is true, an exception will be
thrown if the parameter file is not given.  If this behavior is not
wanted, the method can be overriden.

To check for the legality of the command line, in simple cases it will
suffice to mark some parameters as `required()`.  In more complex
cases a [[parse_command_line]] can be written that calls this one to
parse the command line and then checks that all required command-line
parameters have been read and are consistent.

If a usage error is detected, show_usage() is called and a Usage_error
exception is thrown.

\param[in] argc,argv      Argument count and values as passed to main()
*/
void SimulationCL::parse_command_line(int argc,char *argv[])
{
  std::string parameter_file;

  try {

    CLParameters::parse_command_line(argc,argv,true);

    if (value("help").as<bool>()) {
      show_usage();
      throw Early_stop();
    }
    if (value("version").as<bool>()) {
      show_version();
      throw Early_stop();
    }
    if (value("list-parameters").as<bool>()) {
      show_parameters(std::cerr);
      throw Early_stop();
    }
    if (count("parameter_file")>0) {
	parameter_file=value("parameter_file").as<std::string>();
	parse(parameter_file.c_str());
    } else {
      notify();
    }

  } catch (po::too_many_positional_options_error& e) {
    throw Usage_error();
  } catch (po::invalid_command_line_syntax& e) {
    throw Usage_error();
  } catch (po::invalid_command_line_style& e) {
    throw Usage_error();
  } catch (po::required_option& e) {
    throw Required_parameter_absent(e.what());
  } catch (po::reading_file& e) {
    throw Runtime_error("Cannot read parameter file "+parameter_file);
  }
}

/**
This is automatically called bye parse_command_line() on detecting a
the `--help` option.  You can call it explicitly on catching a
Usage_error excpetion.  Override if the help messega does not apply.
*/
void SimulationCL::show_usage()
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
    << "  --version               Show this version\n"
    << "  --list-parameters       Show accepted parameters\n"
	 ;
}

void SimulationCL::show_version()
{
  std::cerr << name_and_ver << " (with glsim " << VERSION << ")\n";
  std::cerr << copyr << '\n';
  // std::cerr << "glsim is Copyright (C) 2009-2015 Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar\n";
}

/*****************************************************************************
 *
 * class UtilityCL
 *
 */

/**
The constructor defines the common options `--help`, `--version`,
`--verbose`, and `--terse`, but no positional options.

\param[in] name The name of the utility for show_version().
 */
UtilityCL::UtilityCL(const char* name,const char *scope) :
  CLParameters(scope),
  name_and_ver(name)

{
  command_line_options().add_options()
    ("help,h",po::bool_switch(),"help with usage")
    ("verbose,V",po::bool_switch(),"be verbose")
    ("terse,T",po::bool_switch(),"be terse: print as little as possible (e.g. omit headers)")
    ("version",po::bool_switch(),"print version and exit")
    ;
}

/**
Parses the command line (through CLParameters::parse_comand_line())
and act upon the help options (`--help` and `--version`),
displaying the requested help and throwing an
exception to request an early stop.

To check for the legality of the command line, in simple cases it will
suffice to mark some parameters as `required()`.  In more complex
cases a parse_command_line can be written that calls this one to
parse the command line and then checks that all required command-line
parameters have been read and are consistent.

If a usage error is detected, show_usage() is called and a Usage_error
exception is thrown.  This calls notify(), so that if a parameter file
needs to be read it must be placed in another scope, or this function
overridden.

\param[in] argc,argv      Argument count and values as passed to main()
*/
void UtilityCL::parse_command_line(int argc,char *argv[])
{
  try {

    CLParameters::parse_command_line(argc,argv,true);

    if (value("help").as<bool>()) {
      show_usage();
      throw Early_stop();
    }
    if (value("version").as<bool>()) {
      show_version();
      throw Early_stop();
    }
    notify();

  } catch (po::too_many_positional_options_error& e) {
    throw Usage_error();
  } catch (po::invalid_command_line_syntax& e) {
    throw Usage_error();
  } catch (po::invalid_command_line_style& e) {
    throw Usage_error();
  } catch (po::required_option& e) {
    throw Usage_error();
  }
}

void UtilityCL::show_version()
{
  std::cerr << name_and_ver << " (part of glsim " << VERSION << ")\n";
  std::cerr << "Copyright (C) 2009-2015 Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>\n";
}

int UtilityEC(int argc,char *argv[],void (*wmain)(int,char**))
{
  enum return_codes
  {no_error=0,early_stop=1,usage_error=2,runtime_error=10,
   logic_error=20,other_exception=255} ;

  return_codes rcode=no_error;

  try {

    wmain(argc,argv);
  
  } catch (const glsim::Early_stop &e) {
    rcode=early_stop;
  } catch (const glsim::Usage_error &e) {
    std::cerr << "Wrong usage, try --help\n";
    rcode=usage_error;
  } catch (const glsim::Runtime_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
    rcode=runtime_error;
  } catch (const glsim::Logic_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
    rcode=logic_error;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << '\n';
    rcode=other_exception;
  }
  return rcode;
}


} /* namespace */

