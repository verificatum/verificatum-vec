
# Copyright 2008-2018 Douglas Wikstrom
#
# This file is part of Verificatum Elliptic Curve library (VEC).
#
# VEC is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# VEC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General
# Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with VEC. If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.63])
AC_INIT([verificatum-vec],
        [m4_esyscmd([m4 .version.m4])],
        [info@verificatum.com])
AC_LANG_PUSH([C])
AM_INIT_AUTOMAKE([foreign tar-ustar -Wall -Werror -Wno-portability])
AC_CONFIG_SRCDIR([array_alloc_init.c])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

# Disable any default flags set by AC_PROG_CC. This still allows the
# user to specify CFLAGS.
CFLAGSTEST=`echo ${CFLAGS} | tr -d " "`
if test x${CFLAGSTEST} = x;
then
   # This is needed since AC_PROG_CC is implemented incorrectly. It
   # checks if CFLAGS is undefined or defined to the empty string. It
   # should check if CFLAGS is white space.
   CFLAGS="-fmessage-length=0"
fi 

# Checks for programs.
AC_PROG_CC

# We are building a library.
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB(gmp, __gmpz_init, ,
       [AC_MSG_ERROR(["GNU MP library not found, see http://gmplib.org/"])])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])
AC_CHECK_HEADERS([gmp.h], ,
       [AC_MSG_ERROR(["GNU MP header not found, see http://gmplib.org/"])])

gcc extract_GMP_CFLAGS.c -o extract_GMP_CFLAGS

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile])
AC_OUTPUT



