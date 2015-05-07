/*
 * graph_test.cc -- Driver for the lattice classes tests
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

#include <iostream>

#include "test_exception.hh"

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARD_ERROR 99

extern void run_tests();

int main(int argc, char *argv[])
{
  int result=TEST_SUCCESS;
    
  try {
    run_tests();
  } catch(const Test_failure &e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
    result=TEST_FAILURE;
  } catch(const Test_skip &e) {
    std::cerr << e.what() << '\n';
    result=TEST_SKIPPED;
  } catch (const glsim::Runtime_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
    result=TEST_HARD_ERROR;
  } catch (const glsim::Logic_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
    result=TEST_HARD_ERROR;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << '\n';
    result=TEST_HARD_ERROR;
  }

  return result;
}
