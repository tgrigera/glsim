/*
 * gs_1stneighbours.cc -- Compute distribution of first neighbours
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

/** \file gs_1stneighbours.cc
    \ingroup Structure

    \brief Compute the distribution of first-neighbour distance

*/

#include "log.hh"
#include "parameters.hh"
#include "olconfiguration.hh"
#include "avevar.hh"
#include "histogram.hh"
#include "nneighbours.hh"

struct optlst {
public:
  std::vector<std::string> ifiles;
  int    Nbins;
  bool   mt;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("gs_gr")
{
  hidden_command_line_options().add_options()
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),
     "input files")
    ;
  command_line_options().add_options()
    ("Nbins,N",po::value<int>(&options.Nbins)->required(),
     "Number of bins to use")
    ("multi-threaded,m",po::bool_switch(&options.mt)->default_value(false),
     "use multithreaded algorithm for neighbour iteration")
    ;

  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options]  ifile [ifile ....]\n\n"
    << "Computes the distribution of first (nearest) neighbours over the\n"
    << "given trajectories.  You must give -N.\n"
    << "\n"
    << " Options:\n";
  show_command_line_options(std::cerr);
}

/*****************************************************************************/

void process_conf(glsim::OLconfiguration &conf,glsim::AveVar<false> &av,
		  glsim::Histogram& nnd)
{
  static glsim::Subcells NN(nnd.max()/5);
  double d,mind;

  NN.rebuild(conf);
  for (int i=0; i<conf.N; ++i) {
    mind=2*nnd.max()*nnd.max();
    for (auto n=NN.neighbours_begin(i); n!=NN.neighbours_end(i); ++n)
      if ( (d=conf.distancesq(i,*n))<mind) mind=d;
    mind=sqrt(mind);
    nnd.push(mind);
    av.push(mind);
  }
}

void process_conf_mt(glsim::OLconfiguration &conf,glsim::AveVar<false> &av,
		  glsim::Histogram& nnd)
{
  static glsim::Subcells NN(nnd.max()/5);

  NN.rebuild(conf);
  #pragma omp parallel for schedule(static)
  for (int i=0; i<conf.N; ++i) {
    double d,mind=2*nnd.max()*nnd.max();
    for (auto n=NN.neighbours_begin(i); n!=NN.neighbours_end(i); ++n)
      if ( (d=conf.distancesq(i,*n))<mind) mind=d;
    mind=sqrt(mind);
    #pragma omp critical
    {
      nnd.push(mind);
      av.push(mind);
    }
  }
}



/*****************************************************************************
 *
 * main
 *
 */

void wmain(int argc,char *argv[])
{
  glsim::logs.set_stream(std::cout,glsim::error);

  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  void (*process_conf_f)(glsim::OLconfiguration &conf,glsim::AveVar<false> &av,
			 glsim::Histogram& nnd);

  if (options.mt) process_conf_f=process_conf_mt;
  else process_conf_f=process_conf;


  ifs.read();
  double diag=sqrt(conf.box_length[0]*conf.box_length[0] +
	    conf.box_length[1]*conf.box_length[1] +
	    conf.box_length[2]*conf.box_length[2]);
  glsim::Histogram nnd(options.Nbins,0.,diag/2);
  nnd.probability_output();
  glsim::AveVar<false> av;

  do {
    process_conf_f(conf,av,nnd);
  } while (ifs.read());
  std::cout << "# Average nearest neighbour distance = " << av.ave() << '\n';
  std::cout << "# Average nearest neighbour standard dev = " << sqrt(av.var()) << '\n';
  std::cout << "#\n";
  if (nnd.outliers()>0)
    std::cout << "#\n# WARNING: There were " << nnd.outliers() << " outliers in histogram\n#\n";
  std::cout << nnd;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
