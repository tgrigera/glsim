/*
 * bethe_test.hh -- testing for classes defined in bethe.hh
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

#include "bethe.hh"
#include "test_exception.hh"

#include "graph_test.hh"

class BetheLattice_test : public Graph_test<glsim::BetheLattice<double>> {
 public:
  BetheLattice_test(int z,int L);

  void test();
  void test_access();
  void test_iterator();

 private:
  int    z,L;
  glsim::BetheLattice<double> lattice;
} ;

BetheLattice_test::BetheLattice_test(int z_,int L_) :
  z(z_),
  L(L_),
  lattice(z,L)
{}

void BetheLattice_test::test_access()
{
  for (int i=0; i<L; ++i)
    lattice.data()[i]=55.5;

  check_result("test_access",lattice[L-1],55.5);
}

void BetheLattice_test::test_iterator()
{
  // The bidirectional iterator operator are tested by the base class
  Graph_test::test_iterator(lattice);

  // now we test the neighbour stuff
  glsim::BetheLattice<double>::node_iterator ni(lattice,&(lattice[3]));
  *ni=0;
  ni.neighbour(0)=1;
  ni.neighbour(1)=2;
  ni.neighbour(2)=3;
  ni.neighbour(3)=4;

  struct ntest {
    void operator()(double &s)
    {ok = ok && s==n; n+=1.;}
    bool ok;
    double n;

    ntest() : ok(true), n(1.) {}
  } ;

  check_result("node iterator neighbour access",
	       for_each_neighbour(ni,ntest()).ok,true);

  ni.to_neighbour(1);
  check_result("node iterator neighbour movement",*ni,2.);
  check_result("node iterator neighbour movement",ni.neighbour(0),0.);
}

void BetheLattice_test::test()
{
  long zm1toL=1;
  for (int i=1; i<=L; i++) zm1toL*=z-1;
  long N = 1 + z*(zm1toL-1)/(z-2);
  test_info(lattice,N,-1);
  test_plain_access(lattice);
  test_access();
  test_iterator();
}

void run_tests()
{
  BetheLattice_test pt(4,5);
  pt.test();
}
