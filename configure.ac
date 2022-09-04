
# Copyright 2008-2019 Douglas Wikstrom
#
# This file is part of Verificatum Elliptic Curve library (VEC).
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use, copy,
# modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

AC_PREREQ([2.69])
AC_INIT([verificatum-vec],[m4_esyscmd([m4 .version.m4])],[info@verificatum.com])
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



