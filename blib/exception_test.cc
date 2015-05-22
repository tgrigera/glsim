/*
 * exception_test.cc -- test/example for exception
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

/** \file exception_test.cc
    \ingroup Test

Here is an example showing how to print a backtrace and source context
information.  To see source lines and function names for the
backtrace, use the `readbt.sh` script.
*/

#include "scontext.hh"
#include "exception.hh"

glsim::Backtrace f2(int a)
{
  glsim::Backtrace bt;
  double f=44*1322*a;
  return bt;
}

glsim::Backtrace f1(int a)
{
  return f2(a);
}

void throws_exceptions()
{
  throw glsim::Unimplemented("ugly feature",HERE);
  // throw glsim::Invalid_operation("reading",HERE);
  errno = 25;
  throw glsim::Clib_error(HERE);
}

int main(int argc, char *argv[])
{

  std::cout << "\n\n\n*** Backtrace test\n\n";
  glsim::Backtrace bt=f1(25);
  std::cout << bt;

  std::cout << "\n\n\n*** Source context test\n\n";
  glsim::Source_context test(HERE),empty;
  std::cout << HERE << "testing context\n";
  std::cout << test << "this is an earlier context\n";
  std::cout << empty << "an empty context\n";
  std::cout << "Backtrace from a Source_context object:\n";
  std::cout << test.backtrace();


  std::cout << "\n\n\n*** Exception test\n\n";
  try {
    throws_exceptions();
  } catch (glsim::Runtime_error& e) {
    std::cerr << "Caught Runtime_error\n";
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
  } catch (glsim::Logic_error& e) {
    std::cerr << "Caught Logic_error\n";
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
  }
  std::cout << "\nContinuing after exceptions.\n\n";

  return 0;
}
