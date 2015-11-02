/*
 * neighbour_test.cc -- tests for structures to find nearest neighbours
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

#include <algorithm>
#include <boost/timer/timer.hpp>

#include "olconfiguration.hh"
#include "random.hh"
#include "nneighbours.hh"
#include "test_exception.hh"

void test_metric_naive(glsim::OLconfiguration &conf)
{
  boost::timer::auto_cpu_timer t;
  
  std::cout << "Testing NeighbourList_naive...";
  std::cout.flush();

  glsim::NeighbourList_naive TNN(conf.box_length[0]/4);
  TNN.rebuild(conf);
  for (int i=0; i<conf.N; ++i) {
    for (auto pj=TNN.neighbours_begin(i); pj!=TNN.neighbours_end(i); ++pj)
      if (std::find(TNN.neighbours_begin(*pj),TNN.neighbours_end(*pj),i) == TNN.neighbours_end(*pj))
	throw Test_failure(std::to_string(i)+" and "+std::to_string(*pj)+
			   " to be mutual neighbours","they are not");
  }
  std::cout << "OK\n";
}

void test_topological_naive(glsim::OLconfiguration &conf)
{
  glsim::TopologicalNeighbours_naive TNN(8);

  std::array<int,3> nnn{8,20,50};
  std::array<int,4> non;
  for (auto n : nnn) {
    boost::timer::auto_cpu_timer t;
    std::cout << "\nTopologicalNearestneighbours_naive with " << conf.N << " particles and " <<
      n << " neighbours\n";
    TNN.rebuild(conf,n);
    non[0]=8;
    non[1]=*(TNN.neighbours_begin(8));
    non[2]=20;
    non[3]=*(TNN.neighbours_begin(20));
    for (auto i : non) {
      std::cout << "     Neighbours for particle " << i << ":\n     ";
      for (auto j=TNN.neighbours_begin(i); j!=TNN.neighbours_end(i); ++j)
	std::cout << *j << " ";
      std::cout << '\n';
    }
    std::cout << '\n';
  }
}

/*****************************************************************************/

void create_random(glsim::OLconfiguration &conf,int N,double boxl)
{
  conf.N=N;
  conf.step=0;
  conf.time=0;
  conf.box_angles[0]=conf.box_angles[1]=conf.box_angles[2]=90.;
  conf.box_length[0]=conf.box_length[1]=conf.box_length[2]=boxl;

  conf.r=new double[conf.N][3];
  glsim::Uniform_real ranx(0,conf.box_length[0]);
  glsim::Uniform_real rany(0,conf.box_length[1]);
  glsim::Uniform_real ranz(0,conf.box_length[2]);
  for (int i=0; i<conf.N; ++i) {
    conf.r[i][0]=ranx();
    conf.r[i][1]=rany();
    conf.r[i][2]=ranz();
  }
}

void run_tests()
{
  glsim::OLconfiguration conf;
  glsim::Random_number_generator RNG(glsim::gsl_rng_mt19937,382930);
  
  create_random(conf,5000,10);

  test_metric_naive(conf);
  test_topological_naive(conf);
}

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARD_ERROR 99

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
