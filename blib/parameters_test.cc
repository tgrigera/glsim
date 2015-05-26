/*
 * parameters_test.cc -- test/example for Parameters
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

/**
   This short program tests and demonstrates how to define new
   parameters to be read from a parameter file or from the
   command-line.
*/

#include "parameters.hh"

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "log.hh"
#include "parameters.hh"

using namespace glsim;

// We now declare parameters intended to be read from a .ini file.
// As shown, many independent Parameters objects can be declared, and
// all parameters will be consolidated and read from the same file.

class filepars : public Parameters {
public:
  filepars();
} ;

filepars::filepars() : 
  Parameters()
{
  parameter_file_options().add_options()
    ("steps",po::value<int>(),"Number of steps")
    ("alpha",po::value<double>()->default_value(1.2),"Weird tuning parameter")
    ;
}

// This is the second class, it also demonstrate the use of sections
// and allows the user to choose scope.

class morepars : public Parameters {
public:
  morepars(const char* scope=Parameters::default_scope);
} ;

morepars::morepars(const char*scope) : Parameters(scope)
{
  parameter_file_options().add_options()
    ("title",po::value<std::string>()->default_value("[no title]"),"Simulation title")
    ("special.steps",po::value<int>(),"Number of special steps")
    ("special.T",po::value<double>()->required(),"Special temperature")
    ("extra_parameters_file",po::value<std::string>(),"Parameters for second scope")
    ;
}

// A third and final class, as we'll see below, we'll use it to place
// parameters in a different scope.  These will be complete
// independent from the above, and will need a second .ini file and
// separate parsing.

class myscopepars : public Parameters {
public:
  myscopepars(const char* scope);
} ;

myscopepars::myscopepars(const char* scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("other.title",po::value<std::string>()->default_value("[other scope I hope]"),"other title")
    ("other.beta",po::value<double>()->required(),"Another weird param")
    ;
}

int main(int argc, char *argv[])
{
  int rcode=0;

  filepars     fpar;
  morepars     mpar;
  SimulationCL cl("Parameter test 1.0","(C) TSG");

  try {

    cl.parse_command_line(argc,argv);

    std::cout << "Parameters read:\n"
	      << "parameter_file = " << cl.value("parameter_file").as<std::string>() << '\n'
	      << "extra_parameter_file = " << cl.value("extra_parameters_file").as<std::string>() << '\n'
	      << "output_file = "
	      << ( cl.count("output_file")>0 ? cl.value("output_file").as<std::string>() : "N/A" ) << '\n'
	      << "steps = " << fpar.value("steps").as<int>() << '\n'
	      << "alpha = " << fpar.value("alpha").as<double>() << '\n'
	      << "title = " << mpar.value("title").as<std::string>() << '\n'
	      << "special.steps = " << mpar.value("special.steps").as<int>() << '\n'
	      << "special.T = " << mpar.value("special.T").as<double>() << '\n'
      ;


    // Let's load the second scope

    morepars secpar("second");
    myscopepars ssecpar("second");

    // the following would generate an unkown parameter exception, showing
    // that scopes are separate:
    //
    // secpar.parse(cl.value("parameter_file").as<std::string>().c_str());
    secpar.parse(cl.value("extra_parameters_file").as<std::string>().c_str());


    std::cout << "\n\nSecond parfile:\n"
	      << "title = " << secpar.value("title").as<std::string>() << '\n'
	      << "special.steps = ";
    if (secpar.count("special.steps") > 0)
      std::cout << secpar.value("special.steps").as<int>();
    else
      std::cout << "undefined";
    std::cout << '\n'
	      << "special.T = " << secpar.value("special.T").as<double>() << '\n'
	      << "other.title = " << ssecpar.value("other.title").as<std::string>() << '\n'
	      << "other.beta = " << ssecpar.value("other.beta").as<double>() << '\n';


  } catch (glsim::Early_stop &e) {
  } catch (glsim::Usage_error &e) {
    cl.show_usage();
    rcode=1;
  } catch (glsim::Runtime_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
    rcode=1;
  } catch (glsim::Logic_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
    rcode=1;
  } catch (std::exception &e) {
    std::cerr << "Exception: " << e.what() << '\n';
    rcode=1;
  }

  return rcode;
}
