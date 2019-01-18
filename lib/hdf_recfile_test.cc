/*
 * hdf_recfile_test.cc -- test of class HDF_record_file
 *
 * This file is part of olglsim, a numerical simulation class library
 * and helper programs.
 *
 * olglsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * olglsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use olglsim to produced published work, or if you redistribute a
 * modified version of olglsim, or code based on olglsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * olglsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding.
 *
 * olglsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

#include "hdf_file.hh"

struct header {
  char         name[200];
  double       box[3];
  float        a;
  int          i;

  header() :
    box{1.,2.,3.}, a(110.), i(25) {strcpy(name,"Test title");}
} head;

struct record {
  double time;
  int    id[10];

  record() : time(1e-1), id{0,1,2,3,4,5,6,7,8,9} {}
} rec;

class RFtest : public glsim::HDF_record_file {
public:
  RFtest(bool old=false);
  void open(const char *fname)
  {HDF_record_file::open(fname,glsim::HDF_record_file::f_replace,"HDF_record_file test");}

  void open_ro(const char *fname)
  {HDF_record_file::open(fname,glsim::HDF_record_file::f_readonly,"HDF_record_file test");}

  void open(const char *fname,const mode& m, const char *tit,
			  version_require vr=version_any,long minver=0)
  {HDF_record_file::open(fname,m,tit,vr,minver);}

  header hbuffer;
  record rbuffer;

private:
  void declare_header_fields(mode);
  void declare_record_fields(mode);
} ;

RFtest::RFtest(bool old) :
  HDF_record_file(old ? 1 : 3)
{
}

void RFtest::declare_header_fields(mode m)
{
  declare_field(f_header,"name",hbuffer.name,200);
  declare_field(f_header,"box",hbuffer.box,3);
  declare_field(f_header,"a radius",&hbuffer.a);
  declare_field(f_header,"run identity",&hbuffer.i);
}

void RFtest::declare_record_fields(mode m)
{
  declare_field(f_record,"time",&rbuffer.time);
  declare_field(f_record,"id",rbuffer.id,10);
}

void create_file()
{
  RFtest file;
  file.open("HDFrectest.dat");
  file.write_header();
  file.write_record(0);
  file.rbuffer.time=1.;
  file.write_record(2);
  file.rbuffer.time=3.;
  file.append_record();
}

void create_file_old_version()
{
  RFtest file(true);
  file.open("HDFrectest_old.dat");
}

void read_file()
{
  RFtest file;
  file.open_ro("HDFrectest.dat");
  file.hbuffer.a=file.hbuffer.i=0;
  file.read_header();
  std::cout << "Read header:\n"
	    << "     name = " << file.hbuffer.name << '\n'
	    << "     box  = " << file.hbuffer.box[0] << ", "
	    << file.hbuffer.box[1] << ", " << file.hbuffer.box[2] << '\n'
	    << "     a    = " << file.hbuffer.a << '\n'
	    << "     i    = " << file.hbuffer.i << "\n\n";

  std::cout << "size = " << file.size() << '\n';

  file.rbuffer.time=0;
  file.rbuffer.id[0]=file.rbuffer.id[9]=0;
  file.read_record(2);

  std::cout << "time " << file.rbuffer.time << "\n";
  std::cout << "id[0], id[9] " << file.rbuffer.id[0] << ','
	    << file.rbuffer.id[9] << "\n";
}

void read_file_old_version()
{
  RFtest file;
  file.open("HDFrectest_old.dat",glsim::HDF_record_file::f_readonly,
	    0,glsim::HDF_record_file::version_min,2);
  // file.open("HDFrectest_old.dat",glsim::HDF_record_file::f_readonly,
  // 	    0,glsim::HDF_record_file::version_exact);
}

int main(int argc, char *argv[])
{
  try {
    
    create_file();
    // create_file_old_version();
    read_file();
    // read_file_old_version();

  } catch (const glsim::Runtime_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
  } catch (const glsim::Logic_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << '\n';
  }

  return 0;
}
