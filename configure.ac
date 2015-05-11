#
# configure.ac -- autoconf script for glsim
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
# glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015
# by Tomas S. Grigera.
# 
#

###############################################################################
#
# Initialization
#
AC_PREREQ([2.67])
AC_INIT([glsim],[1.0.0],[tgrigera@users.sourceforge.net])
AC_CONFIG_SRCDIR([blib/blib.hh])  # Check just one file as safety check
AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AM_INIT_AUTOMAKE
LT_INIT

# Set version numbers for libtool libraries (note that they are
# different from the package version, see
# https://www.sourceware.org/autobook/autobook/autobook_91.html#SEC91
AC_SUBST(LIBGLSIM_VERSION,[1:0:0])

###############################################################################
#
# Checks for programs to needed to build
#
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LN_S
AC_PROG_INSTALL

#
# Check for programs needed to build documentation.  Failure
# will result only in warning and build will proceed w/o documentation

AC_CHECK_PROG(LATEX,latex,latex,[not found])
AC_CHECK_PROG(BIBTEX,bibtex,bibtex,[not found])
AM_CONDITIONAL([HASLATEX],[test x"$LATEX" != x"not found" -a x"$BIBTEX" != x"not found"])
AC_CHECK_PROG(DOXYGEN,doxygen,doxygen,[not found])
AS_IF([test x"$DOXYGEN" = x"not found"],[AC_MSG_WARN([Doxygen not found --- Reference manual will not be built])])
AM_CONDITIONAL([HAS_DOXYGEN],[test x"$DOXYGEN" != x"not found"])
AM_COND_IF([HAS_DOXYGEN],AC_CONFIG_FILES([doc/doxy.conf]))

###############################################################################
#
# Checks for header files.
#
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h unistd.h])
AC_CHECK_HEADER([gsl/gsl_rng.h],,exit)
AC_CHECK_HEADER([boost/program_options.hpp],,exit)
AC_CHECK_HEADER([boost/timer/timer.hpp],,exit)

###############################################################################
#
# Checks for libraries
#

# GSL
AC_CHECK_LIB([gslcblas], [cblas_dsdot],,exit)
AC_CHECK_LIB([gsl], [gsl_rng_alloc],,exit)

# Boost::system
AC_MSG_CHECKING([for working -lboost_system])
SAVED_LIBS=$LIBS
LIBS="-lboost_system $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/system/error_code.hpp>],
  	           [boost::system::error_code e;])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_LIBBOOST_SYSTEM],[1],[Define if Boost::system is available]),
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without boost::system])])

# Boost::program_options
AC_MSG_CHECKING([for working -lboost_program_options])
SAVED_LIBS=$LIBS
LIBS="-lboost_program_options $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/program_options.hpp>],
    [boost::program_options::options_description opt])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_LIBBOOST_PROGRAM_OPTIONS],[1],[Define if Boost::program_options available]) ,
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without boost::program_options])])

# Boost::serialization
AC_MSG_CHECKING([for working -lboost_serialization])
SAVED_LIBS=$LIBS
LIBS="-lboost_serialization $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/archive/binary_oarchive.hpp>
                    #include <fstream>],
    [ std::ofstream os("tst.dat",std::ios::binary);
      boost::archive::binary_oarchive oa(os);])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_LIBBOOST_SERIALIZATION],[1],[Define if Boost::serialization is available]),
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without boost::serialization])])

# Boost::date_time
AC_MSG_CHECKING([for working -lboost_date_time])
SAVED_LIBS=$LIBS
LIBS="-lboost_date_time $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/date_time/posix_time/posix_time.hpp>],
  	           [boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_LIBBOOST_DATE_TIME],[1],[Define if Boost::date_time is available]),
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without boost::date_time])])

# Boost::timer
AC_MSG_CHECKING([for working -lboost_timer])
SAVED_LIBS=$LIBS
LIBS="-lboost_timer $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/timer/timer.hpp>],
  	           [boost::timer::auto_cpu_timer t;])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_LIBBOOST_TIMER],[1],[Define if Boost::timer is available]),
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without boost::timer])])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STAT
AC_FUNC_STRFTIME

# Generates files
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
	         blib/Makefile
                 lattice/Makefile
                 test/Makefile
		 doc/Makefile
		 doc/version.tex])
#AC_CONFIG_FILES([Makefile
#                 base/Makefile
#                 doc/Makefile])
#AC_CONFIG_FILES([doc/version.tex])
AC_OUTPUT
