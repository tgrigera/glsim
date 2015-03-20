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
#include <boost/timer/timer.hpp>

#include <glsim/cerrors.h>

#include "lattice1D.hh"
#include "lattice2D.hh"

#define BIGSIZE 100000000
#define RTBIGSIZE 10000
// #define BIGSIZE 100
// #define RTBIGSIZE 10

void check_1d()
{
  boost::timer::cpu_timer timer;

  std::cout << "\n\n***** Testing Periodic1DLattice\n\n";

  std::cout << "-- Creating, filling and processing simple array of size "
	    << BIGSIZE << "\n";
  timer.start();
  double *l=new double[BIGSIZE];
  *l=10;
  for (double* d=l+1; d!=l+BIGSIZE; )
    *d=*(d++-1)+1.1;

  for (double* d=l; d!=l+BIGSIZE; d++)
    *d+=sin(*d);

  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Creating, filling and processing Periodic1DLattice through data() pointer\n";
  timer.start();
  Periodic1DLattice<double> lat(BIGSIZE);
  double *d=lat.data();
  *d=10;
  for (; d!=lat.data()+lat.size(); d++)
    *d=*(d-1)+1.1;

  for (d=lat.data(); d!=lat.data()+lat.size(); d++)
    (*d)+=sin(*d);

  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Creating, filling and processing Periodic1DLattice through node iterator\n";
  timer.start();
  Periodic1DLattice<double> lat1(BIGSIZE);
  Periodic1DLattice<double>::node_iterator nd=lat.begin();
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
  Periodic1DLattice<double> lat2(10);
  std::ifstream ifs("graph1d.dat");
  lat2.read(ifs);
  timer.stop();
  std::cout << "   Original: size " << lat.size() << " last: " <<
    lat[lat.size()-1] << '\n';
  std::cout << "   Read    : size " << lat2.size() << " last: " <<
    lat2[lat2.size()-1] << '\n';
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Walking and adding neighbours with neighbour_iterator\n";
  timer.start();
  double sum=0;
  for (d=lat.data(); d!=lat.data()+lat.size(); ++d) {
    Periodic1DLattice<double>::neighbour_iterator n=lat.neighbour_begin(d);
    for (n=lat.neighbour_begin(d); n!=lat.neighbour_end(d); ++n)
      sum+=*n;
  }
  std::cout << "   sum = " << sum << '\n';
  timer.stop();
  std::cout << "  " << boost::timer::format(timer.elapsed());

  std::cout << "-- Walking and adding neighbours with for_each\n";
  timer.start();
  struct s {
    double sum;
    s() : sum(0){}
    void operator()(double *d) {sum+=*d;}
  } summ;
  for (nd=lat.begin(); nd!=lat.end(); ++nd) {
    for_each_neighbour(nd,summ);
  }
  std::cout << "   sum = " << summ.sum << '\n';
  timer.stop();
  std::cout << " " << boost::timer::format(timer.elapsed());
}

// void check_sq()
// {
//   boost::timer::cpu_timer timer;

//   std::cout << "Creating, filling and processing array\n";
//   timer.start();
//   double *l=new double[BIGSIZE];
//   *l=10;
//   for (double* d=l+1; d!=l+BIGSIZE; )
//     *d=*(d++-1)+1.1;

//   for (double* d=l; d!=l+BIGSIZE; d++)
//     *d+=sin(*d);

//   timer.stop();

//   std::cout << boost::timer::format(timer.elapsed());

//   std::cout << "Creating, filling and processing SQ lattice through light iterator\n";
//   timer.start();
//   PeriodicSQLattice<double> lat(RTBIGSIZE,RTBIGSIZE);
//   PeriodicSQLattice<double>::iterator d=lat.begin();
//   *d=10;
//   for (; d!=lat.end(); d++)
//     *d=*(d-1)+1.1;

//   for (d=lat.begin(); d!=lat.end(); d++)
//     (*d)+=sin(*d);

//   timer.stop();
//   std::cout << boost::timer::format(timer.elapsed());

//   std::cout << "Creating, filling and processing 2D lattice through node iterator\n";
//   timer.start();
//   PeriodicSQLattice<double> lat1(RTBIGSIZE,RTBIGSIZE);
//   PeriodicSQLattice<double>::node_iterator nd=lat.node_begin();
//   // *nd=10;
//   // for (; nd!=lat.node_end(); )
//   //   *(nd++)+=1.1;

//   // for (nd=lat.node_begin(); nd!=lat.node_end(); ++nd)
//   //    (*nd)+=sin(*nd);

//   timer.stop();
//   std::cout << boost::timer::format(timer.elapsed());

//   std::cout << "Write and read\n";
//   timer.start();
//   std::ofstream f("graph1d.dat");
//   lat.write(f);
//   f.close();
//   PeriodicSQLattice<double> lat2(2,2);
//   std::ifstream ifs("graph1d.dat");
//   lat2.read(ifs);
//   timer.stop();
//   std::cout << "Original: size " << lat.size() << " last: " <<
//     lat[lat.size()-1] << '\n';
//   std::cout << "Read    : size " << lat2.size() << " last: " <<
//     lat2[lat2.size()-1] << '\n';
//   std::cout << boost::timer::format(timer.elapsed());

//   std::cout << "** Walking through neighbours\n";
//   std::cout << "   Using for + neighbour distance vector\n";

//   int* Ndis=new int[lat.size_y()];
//   int* Sdis=new int[lat.size_y()];
//   int* Wdis=new int[lat.size_x()];
//   int* Edis=new int[lat.size_x()];

//   for (int i=0; i<lat.size_y(); i++) {
//     Ndis[i]=1;
//     Sdis[i]=-1;
//   }
//   Sdis[0]=lat.size_y()-1;
//   Ndis[lat.size_y()-1]=-Sdis[0];
//   for (int i=0; i<lat.size_x(); i++) {
//     Edis[i]=lat.size_y();
//     Wdis[i]=-lat.size_y();
//   }
//   Wdis[0]=lat.size_y()*(lat.size_x()-1);
//   Edis[lat.size_x()-1]=-Wdis[0];

//   timer.start();
//   // for (int i=0; i<lat.size_x(); i++)
//   //   for (int j=0; j<lat.size_y(); j++){
//   //     // std::cout << "i , j "  << i <<','<<j <<
//   //     // 	"  i+W, i+E " << i+Wdis[i]<<','<<i+Edis[i]<<
//   //     // 	"  j+N,j+S " << j+Ndis[j]<<','<<j+Sdis[j] << '\n';
//   //     ptrdiff_t p=lat.id(i,j);
//   //     lat[p]+=lat[p+Wdis[i]]+lat[p+Edis[i]]+lat[p+Ndis[j]]+lat[p+Sdis[j]];
//   //   }
//   for (double* p=&(lat[0]); p!=&(lat[0])+lat.size_x()*lat.size_y(); ++p) {
//     ptrdiff_t i,j;
//     lat.posl(lat.id(p),i,j);
//     (*p)+=*(p+Wdis[i])+*(p+Edis[i])+*(p+Ndis[j])+*(p+Sdis[j]);
//   }
//   timer.stop();
//   std::cout << boost::timer::format(timer.elapsed());
//   std::cout << "   Using neighbour iterator\n";
//   timer.start();
//   struct fu {
//     fu(double *&s) :site(s){}
//     double *&site;
//     void operator()(double* n) {(*site)+=*n;}
//   } ff(d);
//   // for (d=lat.begin(); d!=lat.end(); ++d) {
//   for (PeriodicSQLattice<double>::node_iterator n=lat.node_begin();
//        n!=lat.node_end(); ++n) {
//     for_each_neighbour(lat,n,ff);
//     // for (PeriodicSQLattice<double>::neighbour_iterator n=lat.neighbour_begin(d);
//     //  	 n!=lat.neighbour_end(d); ++n) {
//     //    (*d)+=(*n);
//     // }
//     // Unrolling
//     // PeriodicSQLattice<double>::neighbour_iterator n=lat.neighbour_begin(d);
//     // (*d)+=(*n);
//     // ++n;
//     // (*d)+=(*n);
//     // ++n;
//     // (*d)+=(*n);
//     // ++n;
//     // (*d)+=(*n);
//     // }
//   }
//   timer.stop();
//   std::cout << boost::timer::format(timer.elapsed());

//   std::cout << "** Walking through second-neighbours\n";
//   std::cout << "   Using for + neighbour distance vector\n";
//   timer.start();
//   for (int i=1; i<lat.size_x()-1; i++)
//     for (int j=1; j<lat.size_y()-1; j++) {
//       // std::cout << " i,j= " <<i<<','<<j<< " Wdis[i] "<< Wdis[i] << '\n';
//       // lat(i,j)+=
//       // 	lat(i+Wdis[i]+Wdis[i+Wdis[i]],j)+lat(i+Wdis[i],j+Ndis[j])+lat(i+Wdis[i],j
// 								      // +Sdis[j]);
// 	// lat(i+Wdis[i]+Wdis[i+Wdis[i]],j)+lat(i+Wdis[i],j+Ndis[j])+lat(i+Wdis[i],j+Sdis[j])+
// 	// lat(i+Edis[i]+Edis[i+Edis[i]],j)+lat(i+Edis[i],j+Ndis[j])+lat(i+Edis[i],j+Sdis[j])+
// 	// lat(i,j+Ndis[j]+Ndis[j+Ndis[j]])+lat(i+Edis[i],j+Ndis[j])+lat(i+Wdis[i],j+Ndis[j])+
// 	// lat(i,j+Sdis[j]+Sdis[j+Sdis[j]])+lat(i+Edis[i],j+Sdis[j])+lat(i+Wdis[i],j+Sdis[j]);
//     }
//   timer.stop();
//   std::cout << boost::timer::format(timer.elapsed());
//   std::cout << "   Using neighbour iterator\n";
//   timer.start();
//   for (d=lat.begin(); d!=lat.end(); ++d) {
//     for (PeriodicSQLattice<double>::neighbour_iterator n=lat.neighbour_begin(d);
//      	 n!=lat.neighbour_end(d); ++n) {
//       for (PeriodicSQLattice<double>::neighbour_iterator nn=lat.neighbour_begin(n);
// 	   nn!=lat.neighbour_end(n); ++nn) {
//        (*d)+=(*n);
//       }
//     }
//   }
//   timer.stop();
//   std::cout << boost::timer::format(timer.elapsed());
// }

int main(int argc, char *argv[])
{
  //  glsim::random_number_generator(glsim::gsl_rng_mt19937);

  check_1d();
  // check_sq();
  return 0;
}
