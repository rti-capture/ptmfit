dnl process with autoconf to make configure

AC_INIT([ptmfit], [2.4.0], [vipsip@jiscmail.ac.uk])

# foreign stops complaints about a missing README (we use README.md instead)
# and missing INSTALL (the standard Gnu INSTALL is not very useful)
AM_INIT_AUTOMAKE([foreign])

AM_CONFIG_HEADER(config.h)

AC_PROG_CC
AC_PROG_CXX

AC_CHECK_HEADERS([stdlib.h stdio.h math.h string.h])

PKG_CHECK_MODULES(VIPS, vips-cpp)

AC_SUBST(VIPS_INCLUDES)
AC_SUBST(VIPS_LIBS)

AC_OUTPUT(Makefile)

