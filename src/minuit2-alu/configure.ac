#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([minuit-alu], [6.10], [alusiani@gmail.com])

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config/m4)

AC_CONFIG_SRCDIR([src/MnUserFcn.cxx])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE
# AM_INIT_AUTOMAKE(foreign subdir-objects tar-ustar dist-zip)

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

LT_INIT

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset sqrt])

AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 inc/Makefile
                 inc/Math/Makefile
                 inc/Fit/Makefile
                 inc/Minuit2/Makefile
                 src/Makefile
                 test/Makefile
                 test/MnSim/Makefile
                 test/MnTutorial/Makefile])

AC_OUTPUT
