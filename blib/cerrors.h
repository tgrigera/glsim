/*
 * cerrors.h -- error reporting for C programs
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
