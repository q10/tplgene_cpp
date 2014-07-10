dnl Process this file with autoconf to produce a configure script.
dnl

# The initial information about tplgene++ project
AC_INIT([tplgene++],1.0,[bjma@lbl.gov],[tplgene++])
AC_CONFIG_AUX_DIR([.])
AC_PREREQ([2.69])

# Run script to get version number, stored in GIT_VERSION_FILE
AC_SUBST(GIT_VERSION, esyscmd([./src/GIT-VERSION-GEN]))

#This defines a variable that can be used for the interface versioning
AC_SUBST(VERSION_INFO_CORE, [0:0:0])

AC_ENABLE_SHARED(yes)
AC_ENABLE_STATIC(no)

AC_CANONICAL_HOST

printf '
tplgene++: A tplgene rewrite in C++11

'

# automake initialization
AM_INIT_AUTOMAKE([1.10 -Wall no-define foreign])

# Python parameters
AM_PATH_PYTHON
PY_PREFIX=`$PYTHON -c 'import sys ; print sys.prefix'`
PYTHON_LIBS="-lpython$PYTHON_VERSION"
PYTHON_CFLAGS="-I$PY_PREFIX/include/python$PYTHON_VERSION"
AC_SUBST([PYTHON_LIBS])
AC_SUBST([PYTHON_CFLAGS])

# Use C++11
CXXFLAGS="$CXXFLAGS -std=c++11"

# use the C++ compiler for the following checks
AC_LANG(C++)
AC_LANG(C)
AC_PROG_CXX
AC_PROG_CC_C99
AM_PROG_CC_C_O

dnl add the prefix to the library/include paths
TPLPPDEPS="$TPLPPDEPS_DIR"
AC_SUBST(TPLPPDEPS)

OLD_CFLAGS="$CFLAGS"
OLD_CXXFLAGS="$CXXFLAGS"
OLD_CPPFLAGS="$CPPFLAGS"
OLD_LDFLAGS="$LDFLAGS"

# Set the dir for configure macros
AC_CONFIG_MACRO_DIR([m4])


#
# --with configs
#

# boost
AC_ARG_WITH(boost,
        AC_HELP_STRING(
                [--with-boost=PROGRAM],
                [location of the boost libraries]),
        [boost_prog=$withval])
AC_PATH_PROG(BOOST, boost, $boost_prog)

# Boost Libraries
## check if set with --with-boost
if test "x$BOOST" != "x"; then
   AC_SUBST(BOOST_CPPFLAGS,"-I$BOOST/include")
   AC_SUBST(BOOST_LIBS,"-L$BOOST/lib")
   AC_SUBST(CPPFLAGS,"$CPPFLAGS -I$BOOST/include")
   AC_SUBST(LDFLAGS,"$LDFLAGS -L$BOOST/lib")
   AC_SUBST(BOOST,"$BOOST")
else
   AC_SUBST(BOOST,"$TPLPPDEPS/boost")
fi

## check the BOOST library headers
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([boost/foreach.hpp], [HAVELIB_BOOST=1],
    [AC_MSG_ERROR(You need the Boost libraries.)])
AC_LANG_POP([C++])
AC_SUBST(BOOST_CPPFLAGS)
AC_SUBST(BOOST_LIBS)

AC_SUBST([CFLAGS],["$OLD_CFLAGS"])
AC_SUBST([CXXFLAGS],["$OLD_CXXFLAGS"])
AC_SUBST([CPPFLAGS],["$OLD_CPPFLAGS"])
AC_SUBST([LDFLAGS],["$OLD_LDFLAGS"])


# Check for header files
AC_CHECK_HEADERS([fcntl.h])
AC_CHECK_HEADERS([limits.h])
AC_CHECK_HEADERS([stddef.h])
AC_CHECK_HEADERS([sys/socket.h])
AC_CHECK_HEADERS([sys/time.h])
AC_CHECK_HEADERS([syslog.h])
AC_CHECK_HEADER_STDBOOL

# Check for Functions
AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([getcwd])
AC_CHECK_FUNCS([gethostname])
AC_CHECK_FUNCS([memset])
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([realpath])
AC_CHECK_FUNCS([setenv])
AC_CHECK_FUNCS([sqrt])
AC_CHECK_FUNCS([strcasecmp])
AC_CHECK_FUNCS([strncasecmp])
AC_CHECK_FUNCS([strtol])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MKTIME
AC_FUNC_STRTOD
AC_PROG_RANLIB
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# distribute additional compiler and linker flags
# --> set these variables instead of CFLAGS or LDFLAGS
#AC_SUBST([AM_CFLAGS])
#AC_SUBST([AM_LDFLAGS])


dnl Using Lib tool
AM_PROG_AR
LT_INIT

# Figure out the platform for build
AC_MSG_CHECKING(platform)
IS_DARWIN=false
IS_LINUX=false
IS_UNKNOWN=false
case "$host_os" in
	*darwin*)
		AC_MSG_RESULT(darwin)
		IS_DARWIN=true
		;;
	*linux*)
		AC_MSG_RESULT(linux)
		IS_LINUX=true
		;;
	*)
		AC_MSG_RESULT(unknown)
		IS_UNKNOWN=true
		;;
esac
AM_CONDITIONAL(IS_DARWIN, test "$IS_DARWIN" = true)
AM_CONDITIONAL(IS_LINUX, $IS_LINUX)
AM_CONDITIONAL(IS_UNKNOWN, $IS_UNKNOWN)


# files to generate via autotools (.am or .in source files)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  src/Makefile
])

# generate the final Makefile etc.
AC_OUTPUT
