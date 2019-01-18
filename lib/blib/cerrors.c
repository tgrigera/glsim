/*
 * cerrors.c -- error reporting for C programs
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

#include <execinfo.h>
#define _ERRORS_C
#include "cerrors.h"

static char errbuf[2000];

static void print_backtrace_and_exit(int ecode)
{
  void   *ar[50];
  size_t n;

  n=backtrace(ar,50);
  fprintf(stderr,"===== BACKTRACE START ==========\n");
  backtrace_symbols_fd(ar,n,STDERR_FILENO);
  fprintf(stderr,"===== BACKTRACE END ==========\n");
  exit(ecode);
}

void exit_on_system_error(const char *file,int line,const char *func,
			  const char *s,int ecode)
{
  sprintf(errbuf,"%s:%d (%s): %s",file,line,func,s);
  perror(errbuf);
  print_backtrace_and_exit(ecode);
}

void exit_on_error(const char *file,int line,const char *func,
		   const char *s,int ecode)
{
  fprintf(stderr,"%s:%d (%s): %s\n",file,line,func,s);
  print_backtrace_and_exit(ecode);
}
