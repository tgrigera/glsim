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
#include <iostream>

#include "olconfiguration.hh"
#include "random.hh"
#include "nneighbours.hh"
#include "test_exception.hh"

/*
 * A custom find function that accepts different types for first and
 * last
 *
 */

template<class InputIt1,class InputIt2, class T>
InputIt1 find(InputIt1 first, InputIt2 last, const T& value)
{
    for (; first != last; ++first) {
        if (*first == value) {
            return first;
        }
    }
    return first;
}

/*
 * Tests for metric routines
 *
 */

void test_metric_naive(glsim::OLconfiguration &conf)
{
  std::cout << "Testing NeighbourList_naive...";
  std::cout.flush();

  glsim::NeighbourList_naive TNN(conf.box_length[0]/4);
  TNN.rebuild(conf);
  for (int i=0; i<conf.N; ++i) {
    for (auto pj=TNN.neighbours_begin(i), end=TNN.neighbours_end(i); pj!=end; ++pj)
      if (find(TNN.neighbours_begin(*pj),TNN.neighbours_end(*pj),i) ==
	  TNN.neighbours_end(*pj))
	throw Test_failure(std::to_string(i)+" and "+std::to_string(*pj)+
			   " to be mutual neighbours","they are not");
  }
  std::cout << "OK\n";
}


// print pair in order
std::ostream &operator<<(std::ostream& o,const glsim::Subcells::Pair& p)
{
  if (p.first<=p.second)
    o << '(' << p.first << ',' << p.second << ')';
  else
    o << '(' << p.second << ',' << p.first << ')';
  return o;
}

//
// Test generic NeighbourAlgo against NeighbourList_naive
//
template <typename NeighbourAlgo>
void test_metric(glsim::OLconfiguration &conf,std::string name)
{
  double rc=conf.box_length[0]/4;
  double rcsq=rc*rc;

  std::cout << "Testing " << name << "...";
  std::cout.flush();

  NeighbourAlgo NEW(rc);
  NEW.rebuild(conf);

  std::cout << "building naive pair list...";
  std::cout.flush();
  glsim::NeighbourList_naive TNN(rc,0);
  TNN.rebuild(conf);
  
  std::cout << "comparing...";
  std::cout.flush();

  // Test particle-based
  for (int i=0; i<conf.N; ++i) {
    int np=0;
    for (auto pj=NEW.neighbours_begin(i), end=NEW.neighbours_end(i); pj!=end; ++pj) {
      np++;
      if (conf.distancesq(i,*pj)<=rcsq) {
	
	if (find(TNN.neighbours_begin(i),TNN.neighbours_end(i),*pj) ==
	    TNN.neighbours_end(i))
	  throw Test_failure(std::to_string(i)+" and "+ std::to_string(*pj)+
			     " neighbours in both lists",
			     "missing in naive pair list");
      }
      for (auto pj=TNN.neighbours_begin(i), end=TNN.neighbours_end(i); pj!=end; ++pj) {
	if (find(NEW.neighbours_begin(i),NEW.neighbours_end(i),*pj) ==
	    NEW.neighbours_end(i) )
	  throw Test_failure(std::to_string(i)+" and "+ std::to_string(*pj)+
			     " neighbours in both lists",
			     "missing in "+name);
      }
    }
  }

  // Test pair-based
  double dsq;
  int np=0;
  for (auto pp=NEW.pairs_begin(), end=NEW.pairs_end(); pp!=end; ++pp) {
    dsq=conf.distancesq(pp->first,pp->second);
    if (dsq<=rcsq) {
      glsim::NeighbourList_naive::Pair p1(pp->first,pp->second);
      glsim::NeighbourList_naive::Pair p2(pp->second,pp->first);
      bool foundp1 =
	find(TNN.pairs_begin(),TNN.pairs_end(),p1) != TNN.pairs_end();
      bool foundp2 =
	find(TNN.pairs_begin(),TNN.pairs_end(),p2) != TNN.pairs_end();
      if (!foundp1 && !foundp2) 
	throw Test_failure("pair ("+std::to_string(pp->first)+","
			   +std::to_string(pp->second)+") present in both lists",
			   "missing in naive pair list");
    }
    np++;
  }
  for (auto pp=TNN.pairs_begin(), end=TNN.pairs_end(); pp!=end; ++pp) {
    dsq=conf.distancesq(pp->first,pp->second);
    if (dsq<=rcsq) {
      glsim::NeighbourList_naive::Pair p1(pp->first,pp->second);
      glsim::NeighbourList_naive::Pair p2(pp->second,pp->first);
      bool foundp1 =
	find(NEW.pairs_begin(),NEW.pairs_end(),p1) != NEW.pairs_end();
      bool foundp2 =
	find(NEW.pairs_begin(),NEW.pairs_end(),p2) != NEW.pairs_end();
      if (!foundp1 && !foundp2) 
	  throw Test_failure("pair ("+std::to_string(pp->first)+","
			     +std::to_string(pp->second)+") present in both lists",
			     "missing in new algorithm");
    }
    np++;
    if (np%1000==0) std::cout << " (tested " << np << " candidate pairs)\n";
  }

  std::cout << "OK\n";


}

//
// Test generic NeighbourAlgo with for_each_pair
//
template <typename NeighbourAlgo>
void test_metric_for_each(glsim::OLconfiguration &conf,std::string name)
{
  double rc=conf.box_length[0]/4;
  double rcsq=rc*rc;

  std::cout << "Testing for_each_pair with " << name << "...";
  std::cout.flush();

  NeighbourAlgo NEW(rc);
  NEW.rebuild(conf);

  std::cout << "building naive pair list...";
  std::cout.flush();
  glsim::NeighbourList_naive TNN(rc,0);
  TNN.rebuild(conf);
  
  std::cout << "comparing...";
  std::cout.flush();

  // Test pair-based

  std::vector<std::pair<int,int>> pairs_new;

  for_each_pair(NEW,[&pairs_new](int i,int j,double d){pairs_new.push_back(std::pair<int,int>(i,j));});

  double dsq;
  int np=0;
  for (auto pp=pairs_new.begin(), end=pairs_new.end(); pp!=end; ++pp) {
    dsq=conf.distancesq(pp->first,pp->second);
    glsim::NeighbourList_naive::Pair p1(pp->first,pp->second);
    glsim::NeighbourList_naive::Pair p2(pp->second,pp->first);
    bool foundp1 =
      find(TNN.pairs_begin(),TNN.pairs_end(),p1) != TNN.pairs_end();
    bool foundp2 =
      find(TNN.pairs_begin(),TNN.pairs_end(),p2) != TNN.pairs_end();
    if (!foundp1 && !foundp2) 
      throw Test_failure("pair ("+std::to_string(pp->first)+","
			 +std::to_string(pp->second)+") present in both lists",
			 "missing in naive pair list");
    np++;
    if (np%1000==0) std::cout << " (tested " << np << " candidate pairs)\n";
  }
  for (auto pp=TNN.pairs_begin(), end=TNN.pairs_end(); pp!=end; ++pp) {
    dsq=conf.distancesq(pp->first,pp->second);
    if (dsq<=rcsq) {
      glsim::NeighbourList_naive::Pair p1(pp->first,pp->second);
      glsim::NeighbourList_naive::Pair p2(pp->second,pp->first);
      bool foundp1 =
	find(pairs_new.begin(),pairs_new.end(),p1) != pairs_new.end();
      bool foundp2 =
	find(pairs_new.begin(),pairs_new.end(),p2) != pairs_new.end();
      if (!foundp1 && !foundp2) 
	throw Test_failure("pair ("+std::to_string(pp->first)+","
			   +std::to_string(pp->second)+") present in both lists",
			   "missing in new algorithm");
    }
    np++;
    if (np%1000==0) std::cout << " (tested " << np << " candidate pairs)\n";
  }

  std::cout << "OK\n";
}


//
// Test generic NeighbourAlgo with for_each_pair multithreaded
//
typedef std::vector<std::pair<int,int>> plist_t;

class accum_pairs {
public:
  plist_t plist;
  void operator()(int i,int j,double d) {
    plist.push_back(std::pair<int,int>(i,j));
  }
  void reduce(accum_pairs& a)
  {
    plist.insert(plist.end(),a.plist.begin(),a.plist.end());
  }
} ;
    

template <typename NeighbourAlgo>
void test_metric_for_each_mt(glsim::OLconfiguration &conf,std::string name)
{
  double rc=conf.box_length[0]*4;
  double rcsq=rc*rc;

  std::cout << "Testing for_each_pair (multithreaded) with " << name << "...";
  std::cout.flush();

  NeighbourAlgo NEW(rc);
  NEW.rebuild(conf);

  std::cout << "building naive pair list...";
  std::cout.flush();
  glsim::NeighbourList_naive TNN(rc,0);
  TNN.rebuild(conf);
  
  std::cout << "comparing...";
  std::cout.flush();

  // Test pair-based

  plist_t pairs_new=for_each_pair_mt(NEW,accum_pairs()).plist;
  // for_each_pair_mt(NEW,p
  // 		   [&pairs_new](int i,int j,double d)
  // 		   {
  // 		     pairs_new.push_back(std::pair<int,int>(i,j));
  // 		   }
  // 		   );

  double dsq;
  int np=0;
  for (auto pp=pairs_new.begin(), end=pairs_new.end(); pp!=end; ++pp) {
    dsq=conf.distancesq(pp->first,pp->second);
    glsim::NeighbourList_naive::Pair p1(pp->first,pp->second);
    glsim::NeighbourList_naive::Pair p2(pp->second,pp->first);
    bool foundp1 =
      find(TNN.pairs_begin(),TNN.pairs_end(),p1) != TNN.pairs_end();
    bool foundp2 =
      find(TNN.pairs_begin(),TNN.pairs_end(),p2) != TNN.pairs_end();
    if (!foundp1 && !foundp2) 
      throw Test_failure("pair ("+std::to_string(pp->first)+","
			 +std::to_string(pp->second)+") present in both lists",
			 "missing in naive pair list");
    np++;
    if (np%1000==0) std::cout << " (tested " << np << " candidate pairs)\n";
  }
  for (auto pp=TNN.pairs_begin(), end=TNN.pairs_end(); pp!=end; ++pp) {
    dsq=conf.distancesq(pp->first,pp->second);
    if (dsq<=rcsq) {
      glsim::NeighbourList_naive::Pair p1(pp->first,pp->second);
      glsim::NeighbourList_naive::Pair p2(pp->second,pp->first);
      bool foundp1 =
	find(pairs_new.begin(),pairs_new.end(),p1) != pairs_new.end();
      bool foundp2 =
	find(pairs_new.begin(),pairs_new.end(),p2) != pairs_new.end();
      if (!foundp1 && !foundp2) 
	throw Test_failure("pair ("+std::to_string(pp->first)+","
			   +std::to_string(pp->second)+") present in both lists",
			   "missing in new algorithm");
    }
    np++;
    if (np%1000==0) std::cout << " (tested " << np << " candidate pairs)\n";
  }

  std::cout << "OK\n";
}

/*
 * Tests for topological routines
 *
 */

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

/*****************************************************************************
 *
 * Test driver
 *
 */

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
  
  create_random(conf,300,10);

  test_metric_naive(conf);
  test_metric<glsim::MetricNearestNeighbours>(conf,"all pairs enumeration");
  test_metric<glsim::Subcells>(conf,"subcell algorithm");
  test_metric<glsim::NeighbourList_subcells>(conf,"pair list with subcells");

  test_metric_for_each<glsim::MetricNearestNeighbours>(conf,"all pairs enumeration");
  test_metric_for_each<glsim::NeighbourList_naive>(conf,"naive neighbour list");
  test_metric_for_each<glsim::Subcells>(conf,"subcell algorithm");
  test_metric_for_each<glsim::NeighbourList_subcells>(conf,"pair list with subcells");

  test_metric_for_each_mt<glsim::MetricNearestNeighbours>(conf,"pair list with subcells");

  test_topological_naive(conf);
}

/*****************************************************************************
 *
 * main
 *
 */

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
