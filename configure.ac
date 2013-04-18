# configure.ac for glsim
#
#
# This file is part of glsim, a numerical simulation class library and
# helper programs.
#
AC_PREREQ([2.67])
AC_INIT([glsim],[1.0.0],[tgrigera@users.sourceforge.net])
AC_CONFIG_SRCDIR([blib/blib.nw])
AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AM_INIT_AUTOMAKE
LT_INIT

# Checks for programs.
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_CC
AC_PROG_CXX
LT_INIT
AC_PROG_LN_S
AC_PROG_INSTALL
AC_CHECK_PROG(LATEX,latex,latex,[not found])
AC_CHECK_PROG(BIBTEX,bibtex,bibtex,[not found])
AM_CONDITIONAL([HASLATEX],[test x"$LATEX" != x"not found" -a x"$BIBTEX" != x"not found"])
AC_CHECK_PROG(TANGLE,notangle,notangle,[not found])
AC_CHECK_PROG(WEAVE,noweave,noweave,[not found])
AC_CHECK_PROG(DEF,nodefs,nodefs,[not found])
AC_CHECK_PROG(INDEX,noindex,noindex,[not found])
AC_CHECK_PROG(CPIF,cpif,cpif,[not found])
AS_IF([test x"$TANGLE" = x"not found"],[AC_MSG_WARN([****** noweb not found, no source code editing, see README])])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h unistd.h])
AC_CHECK_HEADER([gsl/gsl_rng.h],,exit)
AC_CHECK_HEADER([boost/program_options.hpp],,exit)

#
# Checks for libraries. NOTE that the c++ libs are not properly
# checked, only their existence is. Consider macro AX_CXX_CHECK_LIB from
# sourceforge (guido draheim). For now I'm trying an approach learned from
# http://nerdland.net/2009/07/detecting-c-libraries-with-autotools/
# with boost
#

# GSL
AC_CHECK_LIB([gslcblas], [cblas_dsdot],,exit)
AC_CHECK_LIB([gsl], [gsl_rng_alloc],,exit)

# Boost::program_options
#AC_CHECK_LIB([boost_program_options],main,,exit)
AC_MSG_CHECKING([for working -lboost_program_options])
SAVED_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS -lboost_program_options"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/program_options.hpp>],
    [boost::program_options::options_description opt])],
  AC_MSG_RESULT(yes) ; [LIBS="-lboost_program_options $LIBS"] ; [HAVE_LIBBOOST_PROGRAM_OPTIONS=1],
  AC_MSG_RESULT(no) ; [AC_MSG_FAILURE([cannot build without boost::program_options])])
LDFLAGS=$SAVED_LDFLAGS


# Netcdf
AC_CHECK_LIB([netcdf], nc_create,,exit)
AC_CHECK_LIB([netcdf_c++], nc_create,,exit)

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
                 test/Makefile
		 doc/Makefile
		 doc/version.tex])
#AC_CONFIG_FILES([Makefile
#                 base/Makefile
#                 doc/Makefile])
#AC_CONFIG_FILES([doc/version.tex])
AC_OUTPUT
