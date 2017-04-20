/*
 * poisson_test.hh -- testing for classes defined in poisson.hh
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

#include "poisson.hh"
#include "test_exception.hh"

#include "graph_test.hh"

class MetricPoisson3D_test : public Graph_test<glsim::MetricPoisson3D<double>> {
 public:
  MetricPoisson3D_test(int L,double box[]);

  void test();
  void test_access();
  void test_iterator();
  void test_bond_access();

 private:
  int    N;
  double cutoff;
  glsim::MetricPoisson3D<double> lattice;
} ;

MetricPoisson3D_test::MetricPoisson3D_test(int N_,double boxl[]) :
  N(N_),
  cutoff(boxl[0]*0.2),
  lattice(N,boxl,boxl[0]*0.2)
{}

void MetricPoisson3D_test::test_access()
{
  for (int i=0; i<N; ++i)
    lattice.data()[i]=55.5;

  check_result("test_access",lattice[N-1],55.5);
}

void MetricPoisson3D_test::test_iterator()
{
  // The bidirectional iterator operators are tested by the base class
  Graph_test::test_iterator(lattice);

  // now we test the neighbour stuff
  glsim::MetricPoisson3D<double>::node_iterator ni(lattice,&lattice[0]);
  ni.to(3);
  *ni=0;
  for (int i=0; i<ni.neighbour_size(); ++i)
    ni.neighbour(i)=i+1;

  struct ntest {
    void operator()(double &s)
    {ok = ok && s==n; n+=1.; }
    bool ok;
    double n;

    ntest() : ok(true), n(1.) {}
  } ;

  check_result("node iterator neighbour access",
	       for_each_neighbour(ni,ntest()).ok,true);

  ni.to(3);
  check_result("node iterator movement",*ni,0.);
  ni.to_neighbour(0);
  check_result("node iterator neighbour movement",*ni,1.);
}

void MetricPoisson3D_test::test_bond_access()
{
  int pairs1=0,pairs2=0;
  std::cout << "Testing pairs, cutoff is " << cutoff << '\n';
  for_each_bond(lattice,[this,&pairs1](double &f,double &s){
      std::cout << "Recieved pair with values " << f << ' ' << s << ", distance "
		<< lattice.oconf().distancesq(lattice.id(&f),lattice.id(&s)) << '\n';
      ++pairs1;
    });
  std::cout << "Received total " << pairs1 << " pairs\n";
  std::for_each(lattice.data(),lattice.data()+lattice.size(),
		[this,&pairs2](double &d){pairs2+=lattice.neighbour_size(lattice.id(&d));});
  pairs2/=2;
  std::cout << "Counted " << pairs2 << " pairs from nieghbour lists\n";
  check_result("Pairs",pairs1,pairs2);
}

void MetricPoisson3D_test::test()
{
  test_info(lattice,N,-1);
  test_plain_access(lattice);
  test_access();
  test_iterator();
  test_bond_access();
}

void run_tests()
{
  glsim::Random_number_generator rng(glsim::gsl_rng_mt19937,3932);
  double box[3]={2.,2.,2.};
  MetricPoisson3D_test pt(50,box);
  pt.test();
}
