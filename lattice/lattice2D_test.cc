/*
 * lattice2D_test.hh -- testing for classes defined in lattice2D.hh
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

#include "lattice2D.hh"
#include "test_exception.hh"

#include "graph_test.hh"

class PeriodicSQLattice_test : public Graph_test<glsim::PeriodicSQLattice<double>> {
 public:
  PeriodicSQLattice_test(int L);

  void test();
  void test_access();
  void test_iterator();

 private:
  int    Lx,Ly;
  glsim::PeriodicSQLattice<double> lattice;
} ;

PeriodicSQLattice_test::PeriodicSQLattice_test(int L) :
  Lx(L), Ly(L+3),
  lattice(Lx,Ly)
{}

void PeriodicSQLattice_test::test_access()
{
  for (int i=0; i<Lx*Ly; ++i)
    lattice.data()[i]=55.5;

  check_result("test_access",lattice(Lx-1,Ly-1),55.5);
}

void PeriodicSQLattice_test::test_iterator()
{
  // The bidirectional iterator operator are tested by the base class
  Graph_test::test_iterator(lattice);

  // now we test the neighbour stuff
  glsim::PeriodicSQLattice<double>::node_iterator ni(lattice,&lattice(0,0));
  lattice(0,0)=0;
  lattice(0,1)=1;  // N
  lattice(1,0)=2;  // E
  lattice(0,lattice.size_y()-1)=3;  // S
  lattice(lattice.size_x()-1,0)=4;  // W

  struct ntest {
    void operator()(double &s)
    {ok = ok && s==n; n+=1.;}
    bool ok;
    double n;

    ntest() : ok(true), n(1.) {}
  } ;

  check_result("node iterator neighbour access",
	       for_each_neighbour(ni,ntest()).ok,true);

  ni.to(lattice.size_x()-1,0);
  check_result("node iterator movement",*ni,4.);
  ni.to_neighbour(1);
  check_result("node iterator neighbour movement",*ni,0.);
}

void PeriodicSQLattice_test::test()
{
  test_info(lattice,Lx*Ly,4);
  test_plain_access(lattice);
  test_access();
  test_iterator();
}

void run_tests()
{
  PeriodicSQLattice_test pt(20);
  pt.test();
}
