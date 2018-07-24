dnl Process this file with autoconf to produce a configure script.
AC_INIT([pmpar], [1.8], [Walter Brisken <wbrisken@lbo.us>])
AM_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([foreign])

dnl Checks for programs.
AC_PROG_CXX

dnl Checks for libraries.
AC_CHECK_LIB(m, erf,,[AC_MSG_ERROR("need libm")])

AC_OUTPUT(Makefile)
