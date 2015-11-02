/*
 * test_exception.hh -- Exception class to report test results
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

#ifndef TEST_EXCEPTION_HH
#define TEST_EXCEPTION_HH

#include "exception.hh"

class Test_failure : public glsim::Logic_error {
public:
  explicit Test_failure(const std::string& expected,const std::string& actual,
			const glsim::Source_context &c=glsim::Source_context() ):
    Logic_error("test FAILED:\n  -- expected: "+expected+
		"\n  --      got: "+actual,c)
  {}
} ;

class Test_skip : public glsim::Logic_error {
public:
  explicit Test_skip(const std::string& why) :
    Logic_error("test SKIPPED: "+why)
  {}
} ;

template <typename resT>
void check_result(const char* testing,resT actual,resT expected)
{
  std::cout << "Testing " << testing << ": ";
  if (actual==expected) {
    std::cout << "OK\n";
    return;
  }
  std::string exp=std::to_string(expected);
  std::string act=std::to_string(actual);
  throw Test_failure(exp,act,HERE);
}

#endif /* TEST_EXCEPTION_HH */
