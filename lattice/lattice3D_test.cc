/*
 * lattice3D_test.hh -- testing for classes defined in lattice3D.hh
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

#include "lattice3D.hh"
#include "test_exception.hh"

#include "graph_test.hh"

class PeriodicSCLattice_test : public Graph_test<glsim::PeriodicSCLattice<double>> {
 public:
  PeriodicSCLattice_test(int L);

  void test();
  void test_access();
  void test_iterator();

 private:
  int    Lx,Ly,Lz;
  glsim::PeriodicSCLattice<double> lattice;
} ;

PeriodicSCLattice_test::PeriodicSCLattice_test(int L) :
  Lx(L), Ly(L+3), Lz(L+1),
  lattice(Lx,Ly,Lz)
{}

void PeriodicSCLattice_test::test_access()
{
  for (int i=0; i<Lx*Ly*Lz; ++i)
    lattice.data()[i]=55.5;

  check_result("test_access",lattice(Lx-1,Ly-1,Lz-1),55.5);
}

void PeriodicSCLattice_test::test_iterator()
{
  // The bidirectional iterator operator are tested by the base class
  Graph_test::test_iterator(lattice);

  // now we test the neighbour stuff
  glsim::PeriodicSCLattice<double>::node_iterator ni(lattice,&lattice(3,3,0));
  lattice(3,3,0)=0;
  lattice(3,4,0)=1;  // N
  lattice(4,3,0)=2;  // E
  lattice(3,2,0)=3;  // S
  lattice(2,3,0)=4;  // W
  lattice(3,3,1)=5;  // U
  lattice(3,3,lattice.size_z()-1)=6;  // D

  struct ntest {
    void operator()(double &s)
    {ok = ok && s==n; n+=1.; }
    bool ok;
    double n;

    ntest() : ok(true), n(1.) {}
  } ;

  check_result("node iterator neighbour access",
	       for_each_neighbour(ni,ntest()).ok,true);

  ni.to(3,3,lattice.size_z()-1);
  check_result("node iterator movement",*ni,6.);
  ni.to_neighbour(4);
  check_result("node iterator neighbour movement",*ni,0.);
}

void PeriodicSCLattice_test::test()
{
  test_info(lattice,Lx*Ly*Lz,6);
  test_plain_access(lattice);
  test_access();
  test_iterator();
}

void run_tests()
{
  PeriodicSCLattice_test pt(20);
  pt.test();
}
