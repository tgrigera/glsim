/*
 * cerrors.h -- error reporting for C programs
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

#ifndef CERRORS_H
#define CERRORS_H

/** \file cerrors.h Error reporting in C
    \ingroup Error

For the occasional C program we provide a simple way to abort the
program with a backtrace through the `ERROR_EXIT` and `SYSERROR_EXIT`
macros.  Also, the `I_AM_HERE` macro is provided as an easy way to
report source position (to stderr), useful when debugging.

*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>

/** Exit with error code ecode and message s, reporting context */
#define ERROR_EXIT(s,ecode) \
        exit_on_error(__FILE__,__LINE__,__func__,s,ecode)

/** Like ERROR_EXIT but also print library error (with perror()) */
#define SYSERROR_EXIT(s,ecode) \
        exit_on_system_error(__FILE__,__LINE__,__func__,s,ecode)

/** Print source context (function and source line number */
#define I_AM_HERE \
  {fprintf(stderr,"(DD) Reached %s:%d (in function %s)\n",__FILE__,__LINE__, \
	  __func__); fflush(stderr); }

#ifdef _ERRORS_C
#define EXTERN
#else
#define EXTERN extern
#endif

#ifdef __cplusplus
#include <cstdio>
extern "C" {
#endif
EXTERN void exit_on_system_error(const char *file,int line,const char *func,
				 const char *s,int ecode);
EXTERN void exit_on_error(const char *file,int line,const char *func,
			  const char *s,int ecode);
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CERRORS_H */
