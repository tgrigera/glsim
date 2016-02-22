/*
 * olconf_create.ccc -- create various types of off-lattice configurations
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

/** \file olconf_create.cc
    \ingroup Offlattice

  gs_olconf_create can be used to create an off-lattice configuration
  and save it in the format (HDF5) used by OLconfiguration_file.

  BUG: All particles have the same mass --- this must be fixed.

 */

#include <assert.h>

#include "parameters.hh"
#include "random.hh"
#include "cerrors.h"
#include "olconfiguration.hh"

struct scomp {
  int Nt;
  int *N;
  double boxl[3];

  scomp() : N(0) {}
  ~scomp() {if (N) delete[] N;}
} ;

void create_random(glsim::OLconfiguration &conf,scomp &SC)
{
  conf.N=0;
  for (int c=0; c<SC.Nt; ++c) conf.N+=SC.N[c];
  conf.step=0;
  conf.time=0;
  conf.box_angles[0]=conf.box_angles[1]=conf.box_angles[2]=90.;
  memcpy(conf.box_length,SC.boxl,3*sizeof(double));

  int i;
  conf.id=new short[conf.N];
  for (i=0; i<conf.N; conf.id[i]=i++) ;
  if (SC.Nt>1) {
    conf.type = new short[conf.N];
    i=0;
    for (int c=0; c<SC.Nt; ++c)
      for (int j=0; j<SC.N[c]; ++j)
	conf.type[i++]=c;
  }
  conf.r=new double[conf.N][3];
  glsim::Uniform_real ranx(0,conf.box_length[0]);
  glsim::Uniform_real rany(0,conf.box_length[1]);
  glsim::Uniform_real ranz(0,conf.box_length[2]);
  for (i=0; i<conf.N; ++i) {
    conf.r[i][0]=ranx();
    conf.r[i][1]=rany();
    conf.r[i][2]=ranz();
  }
}

void create_given_coordinates(glsim::OLconfiguration &conf,scomp &SC)
{
  conf.N=0;
  for (int c=0; c<SC.Nt; ++c) conf.N+=SC.N[c];
  conf.step=0;
  conf.time=0;
  conf.box_angles[0]=conf.box_angles[1]=conf.box_angles[2]=90.;
  memcpy(conf.box_length,SC.boxl,3*sizeof(double));

  int i;
  conf.id=new short[conf.N];
  conf.type = new short[conf.N];
  for (i=0; i<conf.N; conf.id[i]=i++) ;
  i=0;
  for (int c=0; c<SC.Nt; ++c)
    for (int j=0; j<SC.N[c]; ++j)
      conf.type[i++]=c;
  conf.r=new double[conf.N][3];
  for (i=0; i<conf.N; ++i) {
    std::cin >> conf.r[i][0] >> conf.r[i][1] >> conf.r[i][2];
  }
}

void set_velocities(glsim::OLconfiguration &conf,double kT,double mass)
{
  delete[] conf.v;
  conf.v = new double[conf.N][3];
  glsim::Gaussian_distribution gauss(sqrt(kT/mass),0);
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0]=gauss();
    conf.v[i][1]=gauss();
    conf.v[i][2]=gauss();
  }
  
}


/*****************************************************************************
 *
 * options and main
 *
 */

static struct oops {
  std::string   ofile;
  bool          from_given_coordinates;

  unsigned long seed;
  int           N;
  double        density;
  double        boxx,boxy,boxz;
  double        vtemp;
  double        mass;

  oops() :
    N(0),
    density(-1)
  {}
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : UtilityCL("gs_olconf_create")
{
  command_line_options().add_options()
    ("boxx",po::value<double>(&options.boxx)->default_value(-1.),
     "Box length in the X direction")
    ("boxy",po::value<double>(&options.boxy)->default_value(-1.),
     "Box length in the Y direction, -1 means use --boxx")
    ("boxz",po::value<double>(&options.boxz)->default_value(-1.),
     "Box length in the Z direction; -1 means use --boxx")
    ("density,d",po::value<double>(&options.density),"Average density.  If given --boxx etc are ignored")
    ("coordinates-from-stdin",po::bool_switch(&options.from_given_coordinates)->default_value(false),"Read coordinates from stdin")
    ("seed,S",po::value<unsigned long>(&options.seed)->default_value(0),"Random number seed")
    ("mass,m",po::value<double>(&options.mass)->default_value(1.),"Particle mass")
    ("velocities,v",po::value<double>(&options.vtemp)->default_value(0.),
     "Generate Maxwellian velocities with kT=arg (if not given, velocities are not written)")
    ;
  hidden_command_line_options().add_options()
    ("out_file",po::value<std::string>(&options.ofile)->required(),"output file")
    ("Nparts",po::value<int>(&options.N)->required(),"total number of particles")
    ;
  positional_options().add("Nparts",1).add("out_file",1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << " Nparts outfile\n\n"
    << "Create an off-lattice configuration with Nparts total particles and save to outfile.  Options:\n\n";
  show_command_line_options(std::cerr);
}

void wmain(int argc,char *argv[])
{
  CLoptions opt;
  opt.parse_command_line(argc,argv);

  glsim::OLconfiguration conf;

  scomp SC;
  
  SC.Nt=1;
  SC.N=new int[1];
  SC.N[0]=options.N;
  if (options.density>0) {
    double volume=options.N/options.density;
    SC.boxl[0]=pow(volume,1./3.);
    SC.boxl[2]=SC.boxl[1]=SC.boxl[0];
  } else {
    if (options.boxx<0) {
      std::cerr << "Must give -d or -boxx\n";
      throw glsim::Usage_error();
    }
    SC.boxl[0]=options.boxx;
    SC.boxl[1]= options.boxy>0 ? options.boxy : options.boxx;
    SC.boxl[2]= options.boxz>0 ? options.boxz : options.boxx;
  }

  glsim::Random_number_generator RNG(glsim::gsl_rng_mt19937,options.seed);

  if (options.from_given_coordinates)
    create_given_coordinates(conf,SC);
  else 
    create_random(conf,SC);

  conf.name="Created by olconf_create";
  if (options.vtemp>0)
    set_velocities(conf,options.vtemp,options.mass);

  conf.save(options.ofile);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
