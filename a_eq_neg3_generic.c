
/* Copyright 2008-2019 Douglas Wikstrom
 *
 * This file is part of Verificatum Elliptic Curve library (VEC).
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include "vec.h"

#include "a_eq_neg3_generic_macros.h"

#include "jmul_template.h"
#include "jmulsw_template.h"
#include "jsmul_template.h"
#include "jfmul_template.h"

void
vec_jdbl_a_eq_neg3_generic_inner(vec_scratch_mpz_t scratch,
                                 mpz_t X3, mpz_t Y3, mpz_t Z3,
                                 vec_curve *curve,
                                 mpz_t X1, mpz_t Y1, mpz_t Z1);

/* Naive version of multiplication. Only used during development.
   void
   vec_jmul_a_eq_neg3_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
   vec_curve *curve,
   mpz_t X, mpz_t Y, mpz_t Z,
   mpz_t scalar)
   {
   vec_jmul_a_eq_neg3_generic_inner(RX, RY, RZ, curve, X, Y, Z, scalar);
   }
*/

void
vec_jmulsw_a_eq_neg3_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                             vec_curve *curve,
                             mpz_t X, mpz_t Y, mpz_t Z,
                             mpz_t scalar)
{
  vec_jmulsw_a_eq_neg3_generic_inner(RX, RY, RZ, curve, X, Y, Z, scalar);
}

void
vec_jsmul_a_eq_neg3_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                            vec_curve *curve,
                            mpz_t *X, mpz_t *Y, mpz_t *Z,
                            mpz_t *scalars,
                            size_t len)
{
  vec_jsmul_a_eq_neg3_generic_inner(RX, RY, RZ, curve, X, Y, Z, scalars, len);
}

vec_jfmul_tab_ptr
vec_jfmul_precomp_a_eq_neg3_generic(vec_curve *curve,
                                    mpz_t X, mpz_t Y, mpz_t Z,
                                    size_t len)
{
  vec_jfmul_tab_ptr ptr;

  ptr.generic =
    (vec_jfmul_tab_generic_inner*)
    malloc(sizeof(vec_jfmul_tab_generic_inner));

  vec_jfmul_init_a_eq_neg3_generic_inner(ptr.generic, curve, len);
  vec_jfmul_prcmp_a_eq_neg3_generic_inner(curve, ptr.generic, X, Y, Z);

  return ptr;
}

void
vec_jfmul_a_eq_neg3_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                            vec_curve *curve,
                            vec_jfmul_tab_ptr ptr,
                            mpz_t scalar)
{
  vec_jfmul_cmp_a_eq_neg3_generic_inner(RX, RY, RZ, curve, ptr.generic, scalar);
}

void
vec_jfmul_free_a_eq_neg3_generic(vec_jfmul_tab_ptr ptr)
{
  vec_jfmul_clear_free_a_eq_neg3_generic_inner(ptr.generic);
}

void
vec_jdbl_a_eq_neg3_generic(vec_scratch_mpz_t scratch,
                           mpz_t X3, mpz_t Y3, mpz_t Z3,
                           vec_curve *curve,
                           mpz_t X1, mpz_t Y1, mpz_t Z1)
{
  vec_jdbl_a_eq_neg3_generic_inner(scratch, X3, Y3, Z3, curve, X1, Y1, Z1);
}

void
vec_jadd_a_eq_neg3_generic(vec_scratch_mpz_t scratch,
                           mpz_t X3, mpz_t Y3, mpz_t Z3,
                           vec_curve *curve,
                           mpz_t X1, mpz_t Y1, mpz_t Z1,
                           mpz_t X2, mpz_t Y2, mpz_t Z2)
{
  vec_jadd_generic(scratch, X3, Y3, Z3, curve, X1, Y1, Z1, X2, Y2, Z2);
}


/* long */
/* time_jdbl_a_eq_neg3_generic(int test_time, mpz_t X, mpz_t Y) */
/* { */
/*   long i; */
/*   int t; */
/*   felem x; */
/*   felem y; */
/*   felem z; */

/*   mpz_t_to_felem(x, X); */
/*   mpz_t_to_felem(y, Y); */
/*   felem_one(z); */

/*   t = clock(); */

/*   i = 0; */
/*   do */
/*     { */
/*       point_double(x, y, z, x, y, z); */
/*       i++; */
/*     } */
/*   while (!vec_done(t, test_time)); */

/*   return i; */
/* } */

/* long */
/* time_jadd_a_eq_neg3_generic(int test_time, mpz_t X, mpz_t Y) */
/* { */
/*   long i; */
/*   int t; */

/*   felem rx; */
/*   felem ry; */
/*   felem rz; */

/*   felem x; */
/*   felem y; */
/*   felem z; */

/*   mpz_t_to_felem(x, X); */
/*   mpz_t_to_felem(y, Y); */
/*   felem_one(z); */

/*   point_double(rx, ry, rz, x, y, z); */

/*   t = clock(); */

/*   i = 0; */
/*   do */
/*     { */
/*       point_add(rx, ry, rz, rx, ry, rz, 0, x, y, z); */
/*       i++; */
/*     } */
/*   while (!vec_done(t, test_time)); */

/*   return i; */
/* } */
