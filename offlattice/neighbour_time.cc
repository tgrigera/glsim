/*
 * neighbour_time.cc -- time different nearest-neighbour algorithms
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
 * by Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>
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
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_duration.hpp>
#include <iostream>

#include "log.hh"
#include "parameters.hh"
#include "olconfiguration.hh"
#include "random.hh"
#include "nneighbours.hh"
#include "test_exception.hh"

/*****************************************************************************
 *
 * options
 *
 */

static struct ooptions {
  int    Nparticles;
  double density;
  double rc;

  bool   time_naive,time_list_naive,time_subcell,time_list_subcell;
  bool   multithreaded;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : UtilityCL("gs_neighbour_time")
{
  command_line_options().add_options()
    ("naive",po::bool_switch(&options.time_naive),
     "time naive algorithm (traverse all pairs)")
    ("list_naive",po::bool_switch(&options.time_list_naive),
     "time pair list with naive list building")
    ("subcells",po::bool_switch(&options.time_subcell),
     "time subcell algorithm")
    ("list_subcells",po::bool_switch(&options.time_list_subcell),
     "time pair list with list built through subcell algorithm")
    ("multithreaded,m",po::bool_switch(&options.multithreaded),
     "test also multithreaded version of for_each_pair")
    ;  
  hidden_command_line_options().add_options()
    ("Npart",po::value<int>(&options.Nparticles)->required(),
     "numer of particles")
    ("density",po::value<double>(&options.density)->required(),"density")
    ("rc",po::value<double>(&options.rc)->required(),"cutoff")
    ;
  positional_options().add("Npart",1).add("density",1).add("rc",1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "\nusage: " << progname << " [options] Nparticles density rc\n\n"
    << "Time nearest neighbour algorithms with configurations of given density,\n"
    << "number of particles, and cutoff.  Only algorithms explicitly requested\n"
    << "in the options will be timed.\n"
    << "\nOptions:\n";
  show_command_line_options(std::cerr);
  std::cerr  << "\n";
}

/*****************************************************************************
 *
 * Timing
 *
 */

//
// A generic for_each_pair
//
template <typename Function,typename NeighboursT>
Function for_each_pair_generic(glsim::OLconfiguration& conf,NeighboursT& NN,Function f)
{
  for (auto p = NN.pairs_begin(), end=NN.pairs_end(); p!=end; ++p) {
    double dsq=conf.distancesq(p->first,p->second);
    if (dsq<=NN.cutoffsq()) f(p->first,p->second,dsq);
  }
  return f;
}

struct energy {
  energy() : ener(0) {}
  void operator()(int i,int j,double dist)
  { 
    double ds=dist*dist*dist;
    ener+=ds*ds-2*ds;
  }
  void reduce(energy& e)
  { ener+=e.ener;}
  double ener;
} ;

std::ostream& operator<<(std::ostream& o,const boost::timer::cpu_times& times)
{
  using namespace boost::posix_time;
  using namespace boost::gregorian;
  using namespace boost::timer;
			    
  time_duration total,user,system;

  o << format(times,2,"%t") << " (CPU) \t" << format(times,2,"%w") << " (wall)";

  return o;
}

template <typename NearestT>
double time_algo(NearestT& NN,glsim::OLconfiguration& conf)
{
  boost::timer::cpu_timer     timer,total_timer;
  boost::timer::cpu_times     times;
  double E1=0,E2=0,E3=0;

  total_timer.start();

  timer.start();
  std::cout << "     build structures...  \t"; std::cout.flush();
  NN.rebuild(conf);
  times=timer.elapsed();
  std::cout << times << '\n';
  
  timer.start();
  std::cout << "     loop generic (x10)...\t";  std::cout.flush();
  for (int i=0; i<10; ++i)
    E1+= for_each_pair_generic(conf,NN,energy()).ener;
  times=timer.elapsed();
  std::cout << times << '\n';

  timer.start();
  std::cout << "     loop ad hoc (x10)... \t";   std::cout.flush();
  for (int i=0; i<10; ++i)
    E2 += glsim::for_each_pair(NN,energy()).ener;
  times=timer.elapsed();
  std::cout << times << '\n';

  if (options.multithreaded) {
    timer.start();
    std::cout << "     loop multithread (x10)... \t";   std::cout.flush();
    for (int i=0; i<10; ++i)
      E3 += glsim::for_each_pair_mt(NN,energy()).ener;
      // E3 += for_each_pair_local(NN,energy()).ener;
    times=timer.elapsed();
    std::cout << times << '\n';
  }

  times=total_timer.elapsed();
  std::cout << "     total\t\t\t" << times << '\n';
  std::cout << "     energies = " << E1 << ' ' << E2 << ' ' << E3 << '\n';

  return E1+E2+E3;
}

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

/*****************************************************************************
 *
 * main
 *
 */

void wmain(int argc,char *argv[])
{
  CLoptions opt;
  opt.parse_command_line(argc,argv);

  glsim::OLconfiguration conf;
  glsim::Random_number_generator RNG(glsim::gsl_rng_mt19937,382930);

  glsim::logs.set_stream(std::cerr,glsim::warn);

  double volume=options.Nparticles/options.density;
  double boxl=pow(volume,1./3.);
  create_random(conf,options.Nparticles,boxl);
  double rc=options.rc;

  std::cout << "# All times in seconds\n\n";

  if (options.time_naive) {
    std::cout << "#### Naive (try all pairs)\n";
    glsim::MetricNearestNeighbours NN(rc);
    time_algo(NN,conf);
  }
  if (options.time_list_naive) {
    std::cout << "#### Neighbour list (list built naively)\n";
    glsim::NeighbourList_naive  NLN(rc,rc*0.2);
    time_algo(NLN,conf);
  }
  if (options.time_subcell) {
    std::cout << "#### Subcells\n";
    glsim::Subcells  NS(rc);
    time_algo(NS,conf);
  }
  if (options.time_list_subcell) {
    std::cout << "#### Neighbour list (list build with subcells)\n";
    glsim::NeighbourList_subcells  NLS(rc,rc*0.2);
    time_algo(NLS,conf);
  }
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
