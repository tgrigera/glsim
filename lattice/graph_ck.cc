/*
 * graph_ck.cc -- Checks for the lattice classes
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

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <boost/timer/timer.hpp>

#include <glsim/cerrors.h>

#include <typeinfo>
#include "bethe.hh"
#include "lattice1D.hh"
#include "lattice2D.hh"
#include "lattice3D.hh"


#define BIGSIZE 100000000
#define RTBIGSIZE 10000
#define QRTBIGSIZE 400
// #define BIGSIZE 100
// #define RTBIGSIZE 10

/******************************************************************************/

void check_bethe()
{
  using glsim::BetheLattice;
  
  boost::timer::cpu_timer timer;

  std::cout << "\n\n***** Testing BetheLattice\n\n";

  int z=4,L=14;
  int N=1+z*(powf(z-1,L)-1)/(z-2);

  double *l=new double[N];
  std::cout << "-- Filling and processing simple array of size "
	    << N << "\n";
  timer.start();
  *l=10;
  for (double* d=l+1; d!=l+N; )
    *d=*(d++-1)+1.1;

  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Creating BetheLattice\n";
  timer.start();
  BetheLattice<double> lat(z,L);
  std::cout << "   Bethe size " << lat.size() << '\n';
  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());
  std::cout << "-- Filling and processing BetheLattice through data() pointer\n";
  timer.start();
  double *d=lat.data();
  *d=10;
  for (; d!=lat.data()+lat.size(); d++)
    *d=*(d-1)+1.1;
  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Filling and processing BetheLattice through node iterator\n";
  timer.start();
  BetheLattice<double>::node_iterator nd=lat.begin();
  *nd=10;
  double old=*nd;
  for (; nd!=lat.end(); ++nd) {
    *nd=old+1.1;
    old=*nd;
  }

  for (nd=lat.begin(); nd!=lat.end(); ++nd)
     (*nd)+=sin(*nd);

  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Writing and reading lattice\n";
  timer.start();
  std::ofstream f("graph1d.dat");
  lat.write(f);
  f.close();
  BetheLattice<double> lat2(2,2);
  std::ifstream ifs("graph1d.dat");
  lat2.read(ifs);
  timer.stop();
  std::cout << "   Original: size " << lat.size() << " last: " <<
    lat[lat.size()-1] << '\n';
  std::cout << "   Read    : size " << lat2.size() << " last: " <<
    lat2[lat2.size()-1] << " connectivity " << lat.coordination_number() << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Walking and adding neighbours through iterator\n";
  timer.start();
  double sum=0;
  for (nd=lat.begin(); nd!=lat.end(); ++nd) {
    for (int i=0; i<nd.neighbour_size(); ++i)
      sum+=nd.neighbour(i);
  }
  std::cout << "   sum = " << sum << '\n';
  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Walking and adding neighbours with for_each\n";
  timer.start();
  struct s {
    double sum;
    s() : sum(0){}
    void operator()(double d) {sum+=d;}
  } ;
  double ss=0;
  for (nd=lat.begin(); nd!=lat.end(); ++nd) {
    ss+=for_each_neighbour(nd,s()).sum;
  }
  std::cout << "   sum = " << ss << '\n';
  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());
}


/******************************************************************************/


void check_1d()
{
  using glsim::Periodic1DLattice;
  
  boost::timer::cpu_timer timer;

  std::cout << "\n\n***** Testing Periodic1DLattice\n\n";

  std::cout << "-- Filling and processing simple array of size "
	    << BIGSIZE << "\n";
  double *l=new double[BIGSIZE];
  timer.start();
  *l=10;
  for (double* d=l+1; d!=l+BIGSIZE; ++d)
    *d=*(d-1)+1.1;
  timer.stop();
  std::cout << "   First, second and last: " << l[0] << " " 
	    << l[1] << " " << l[BIGSIZE-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Filling and processing Periodic1DLattice through data() pointer\n";
  Periodic1DLattice<double> lat(BIGSIZE);
  timer.start();
  double *d=lat.data();
  *d=10;
  for (++d; d!=lat.data()+lat.size(); d++)
    *d=*(d-1)+1.1;
  timer.stop();
  std::cout << "   First, second and last: " << lat[0] << " " 
	    << lat[1] << " " << lat[BIGSIZE-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Filling and processing Periodic1DLattice through node iterator\n";
  timer.start();
  Periodic1DLattice<double>::node_iterator nd=lat.begin();
  *nd=10;
  double old=*nd;
  for (++nd; nd!=lat.end(); ++nd) {
    *nd=old+1.1;
    old=*nd;
  }
  timer.stop();
  std::cout << "   First, second and last: " << lat[0] << " " 
	    << lat[1] << " " << lat[BIGSIZE-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());


  Periodic1DLattice<double> lat1(BIGSIZE);


  std::cout << "-- Writing and reading lattice\n";
  timer.start();
  std::ofstream f("graph1d.dat");
  lat.write(f);
  f.close();
  Periodic1DLattice<double> lat2(10);
  std::ifstream ifs("graph1d.dat");
  lat2.read(ifs);
  timer.stop();
  std::cout << "   Original: size " << lat.size() << " last: " <<
    lat[lat.size()-1] << '\n';
  std::cout << "   Read    : size " << lat2.size() << " last: " <<
    lat2[lat2.size()-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Walking and adding neighbours with raw pointers\n";
  timer.start();
  double sum=0;
  d=lat.data();
  sum+=*(d+1)+(*d+lat.size()-1);
  for (; d!=lat.data()+lat.size()-1; ++d)
    sum+=*(d-1)+*(d+1);
  sum+=*(d-1)+*(lat.data());
  std::cout << "   sum = " << sum << '\n';
  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Walking and adding neighbours with node_iterator\n";
  timer.start();
  sum=0;
  for (nd=lat.begin(); nd!=lat.end(); ++nd) {
  //   for (int i=0; i<nd.neighbour_size(); ++i)
  //     sum+=nd.neighbour(i);
    sum+=nd.L()+nd.R();
  }
  std::cout << "   sum = " << sum << '\n';
  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Walking and adding neighbours with for_each\n";
  timer.start();
  struct s {
    double sum;
    s() : sum(0){}
    void operator()(double &d) {sum+=d;}
  } summ;
  for (nd=lat.begin(); nd!=lat.end(); ++nd) {
    // glsim::for_each_neighbour(nd,summ);
  }
  std::cout << "   sum = " << summ.sum << '\n';
  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());
}

void check_sq()
{
  boost::timer::cpu_timer timer;

  std::cout << "\n\n***** Testing PeriodicSQLattice\n\n";

  std::cout << "-- Filling and processing simple array of size "
	    << BIGSIZE << "\n";
  double *l=new double[BIGSIZE];
  timer.start();
  *l=10;
  for (double* d=l+1; d!=l+BIGSIZE; ++d)
    *d=*(d-1)+1.1;
  std::cout << "   First, second and last: " << l[0] << " " 
	    << l[1] << " " << l[BIGSIZE-1] << '\n';

  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Filling and processing PeriodicSQLattice through data() pointer\n";
  glsim::PeriodicSQLattice<double> lat(RTBIGSIZE,RTBIGSIZE);
  timer.start();
  double *d=lat.data();
  *d=10;
  for (++d; d!=lat.end(); ++d)
    *d=*(d-1)+1.1;
  std::cout << "   First, second and last : " << *(lat.data()) << " " 
	    << *(lat.data()+1) << " " << *(--d) << '\n';

  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Filling and processing PeriodicSQLattice through node iterator\n";
  glsim::PeriodicSQLattice<double> lat1(RTBIGSIZE,RTBIGSIZE);
  timer.start();
  glsim::PeriodicSQLattice<double>::node_iterator nd=lat1.begin();
  *nd=10;
  for (++nd; nd!=lat1.end(); ) {
    double& last=*nd;
    ++nd;
    *nd=last +=1.1;
  }
  // for (nd=lat1.begin(); nd!=lat1.end(); ++nd)
  //    (*nd)+=sin(*nd);
  std::cout << "   First, second and last : " << *(lat1.data()) << " " 
	    << *(lat1.data()+1) << " " << *(--nd) << '\n';

  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  {
  std::cout << "-- Write and read\n";
  timer.start();
  std::ofstream f("graph1d.dat");
  lat.write(f);
  f.close();
  std::ifstream ifs("graph1d.dat");
  glsim::PeriodicSQLattice<double> lat2(1,1);
  lat2.read(ifs);
  timer.stop();
  std::cout << "   Original: size " << lat.size() << " last: " <<
    lat[lat.size()-1] << '\n';
  std::cout << "   Read    : size " << lat2.size() << " last: " <<
    lat1[lat1.size()-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());
  }

  std::cout << "-- Walking and adding neighbours using for + neighbour distance vector\n";

  int* Ndis=new int[lat.size_y()];
  int* Sdis=new int[lat.size_y()];
  int* Wdis=new int[lat.size_x()];
  int* Edis=new int[lat.size_x()];
  for (int i=0; i<lat.size_y(); i++) {
    Ndis[i]=1;
    Sdis[i]=-1;
  }
  Sdis[0]=lat.size_y()-1;
  Ndis[lat.size_y()-1]=-Sdis[0];
  for (int i=0; i<lat.size_x(); i++) {
    Edis[i]=lat.size_y();
    Wdis[i]=-lat.size_y();
  }
  Wdis[0]=lat.size_y()*(lat.size_x()-1);
  Edis[lat.size_x()-1]=-Wdis[0];

  lat1=lat;

  timer.start();
  for (double* p=lat.data(); p!=&(lat[0])+lat.size_x()*lat.size_y(); ++p) {
    glsim::id_t i,j;
    glsim::id_t l=lat.id(p);
    i=l/lat.size_y();
    j=l-i*lat.size_y();
    (*p)+=*(p+Wdis[i])+*(p+Edis[i])+*(p+Ndis[j])+*(p+Sdis[j]);
  }
  timer.stop();
  std::cout << "   First, second and last : " << lat[0] << " " 
	    << lat[1] << " " << lat[lat.size()-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Walking and adding neighbours using for_each\n";
  timer.start();
  struct fu {
    fu() : sum(0) {}
    void operator()(double& n) {sum+=n;}
    double sum;
  } fuu;
  for (nd=lat1.begin(); nd!=lat1.end(); ++nd) {
    *nd += for_each_neighbour(nd,fu()).sum;
  }
  timer.stop();
  std::cout << "   First, second and last : " << lat1[0] << " " 
	    << lat1[1] << " " << lat1[lat1.size()-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());
}

///////////////////////////////////////////////////////////////////////////////
//
// Simple cubic

void check_sc()
{
  boost::timer::cpu_timer timer;

  std::cout << "\n\n***** Testing PeriodicSQLattice\n\n";

  int NN=QRTBIGSIZE*QRTBIGSIZE*QRTBIGSIZE;
  std::cout << "-- Filling and processing simple array of size "
	    << NN << "\n";
  double *l=new double[NN];
  timer.start();
  *l=10;
  for (double* d=l+1; d!=l+NN; ++d)
    *d=*(d-1)+1.1;
  std::cout << "   First, second and last: " << l[0] << " " 
	    << l[1] << " " << l[NN-1] << '\n';

  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Filling and processing PeriodicSCLattice through data() pointer\n";
  glsim::PeriodicSCLattice<double> lat(QRTBIGSIZE,QRTBIGSIZE,QRTBIGSIZE);
  timer.start();
  double *d=lat.data();
  *d=10;
  for (++d; d!=lat.end(); ++d)
    *d=*(d-1)+1.1;
  std::cout << "   First, second and last : " << *(lat.data()) << " " 
	    << *(lat.data()+1) << " " << *(--d) << '\n';

  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Filling and processing PeriodicSCLattice through node iterator\n";
  glsim::PeriodicSCLattice<double> lat1(QRTBIGSIZE,QRTBIGSIZE,QRTBIGSIZE);
  timer.start();
  glsim::PeriodicSCLattice<double>::node_iterator nd=lat1.begin();
  *nd=10;
  for (++nd; nd!=lat1.end(); ) {
    double& last=*nd;
    ++nd;
    *nd=last +=1.1;
  }
  std::cout << "   First, second and last : " << *(lat1.data()) << " " 
	    << *(lat1.data()+1) << " " << *(--nd) << '\n';

  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  {
  std::cout << "-- Write and read\n";
  timer.start();
  std::ofstream f("graph1d.dat");
  lat.write(f);
  f.close();
  std::ifstream ifs("graph1d.dat");
  glsim::PeriodicSCLattice<double> lat2(1,1,1);
  lat2.read(ifs);
  timer.stop();
  std::cout << "   Original: size " << lat.size() << " last: " <<
    lat[lat.size()-1] << '\n';
  std::cout << "   Read    : size " << lat2.size() << " last: " <<
    lat1[lat1.size()-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());
  }

  std::cout << "-- Walking and adding neighbours using for + neighbour distance vector\n";

  int* Ndis=new int[lat.size_y()];
  int* Sdis=new int[lat.size_y()];
  int* Wdis=new int[lat.size_x()];
  int* Edis=new int[lat.size_x()];
  int* Udis=new int[lat.size_z()];
  int* Ddis=new int[lat.size_z()];
  for (int i=0; i<lat.size_z(); i++) {
    Udis[i]=1;
    Ddis[i]=-1;
  }
  Ddis[0]=lat.size_z()-1;
  Udis[lat.size_z()-1]=-Ddis[0];
  for (int i=0; i<lat.size_y(); i++) {
    Ndis[i]=lat.size_z();
    Sdis[i]=-lat.size_z();
  }
  Sdis[0]=lat.size_z()*(lat.size_y()-1);
  Ndis[lat.size_y()-1]=-Sdis[0];
  for (int i=0; i<lat.size_x(); i++) {
    Edis[i]=lat.size_y()*lat.size_z();
    Wdis[i]=-lat.size_y()*lat.size_z();
  }
  Wdis[0]=lat.size_z()*lat.size_y()*(lat.size_x()-1);
  Edis[lat.size_x()-1]=-Wdis[0];

  lat1=lat;

  timer.start();
  for (double* p=lat.data(); p!=&(lat[0])+lat.size_x()*lat.size_y()*lat.size_z(); ++p) {
    int i,j,k;
    glsim::id_t l=lat.id(p);
    i=l/(lat.size_y()*lat.size_z());
    j=(l-i*lat.size_y()*lat.size_z())/lat.size_z();
    k=l-i*lat.size_y()*lat.size_z()-j*lat.size_z();
    (*p)+=*(p+Wdis[i])+*(p+Edis[i])+*(p+Ndis[j])+*(p+Sdis[j])+*(p+Udis[k])+
      *(p+Ddis[k]);
  }
  timer.stop();
  std::cout << "   First, second and last : " << lat[0] << " " 
	    << lat[1] << " " << lat[lat.size()-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Walking and adding neighbours using for_each\n";
  timer.start();
  struct fu {
    fu() : sum(0) {}
    void operator()(double& n) {sum+=n;}
    double sum;
  } fuu;
  for (nd=lat1.begin(); nd!=lat1.end(); ++nd) {
    *nd += for_each_neighbour(nd,fu()).sum;
  }
  timer.stop();
  std::cout << "   First, second and last : " << lat1[0] << " " 
	    << lat1[1] << " " << lat1[lat1.size()-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());
}

#define TEST_SQ 2

int main(int argc, char *argv[])
{
  check_bethe();
  check_1d();
  check_sq();
  check_sc();
}
