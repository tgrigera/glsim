/*
 * test_exception.hh -- Exception class to report test results
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use glsim to produced published work, or if you redistribute a
 * modified version of glsim, or code based on glsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * glsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding.  * glsim
 * distribution.
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
