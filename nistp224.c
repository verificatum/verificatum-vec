
/*
 * Copyright 2008-2018 Douglas Wikstrom
 *
 * This file is part of Verificatum Elliptic Curve library (VEC).
 *
 * VEC is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * VEC is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with VEC. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include "vec.h"

#include "ecp_nistp224_core.c"
#include "ecp_nistp224_util.c"

#include "nistp224_macros.h"

#include "jmul_template.h"
#include "jmulsw_template.h"
#include "jsmul_template.h"
#include "jfmul_template.h"

/* Naive version of multiplication. Only used during development.
void
vec_jmul_nistp224(mpz_t RX, mpz_t RY, mpz_t RZ,
                  vec_curve *curve,
                  mpz_t X, mpz_t Y, mpz_t Z,
                  mpz_t scalar)
{
  felem x;
  felem y;
  felem z;
  felem rx;
  felem ry;
  felem rz;

  mpz_t_to_felem(x, X);
  mpz_t_to_felem(y, Y);
  mpz_t_to_felem(z, Z);

  vec_jmul_nistp224_inner(rx, ry, rz,
                          curve,
                          x, y, z,
                          scalar);
  felem_to_mpz_t(RX, rx);
  felem_to_mpz_t(RY, ry);
  felem_to_mpz_t(RZ, rz);
}
*/

void
vec_jmulsw_nistp224(mpz_t RX, mpz_t RY, mpz_t RZ,
                    vec_curve *curve,
                    mpz_t X, mpz_t Y, mpz_t Z,
                    mpz_t scalar)
{
  felem x;
  felem y;
  felem z;
  felem rx;
  felem ry;
  felem rz;

  mpz_t_to_felem(x, X);
  mpz_t_to_felem(y, Y);
  mpz_t_to_felem(z, Z);

  vec_jmulsw_nistp224_inner(rx, ry, rz,
                            curve,
                            x, y, z,
                            scalar);
  felem_to_mpz_t(RX, rx);
  felem_to_mpz_t(RY, ry);
  felem_to_mpz_t(RZ, rz);
}

void
vec_jsmul_nistp224(mpz_t RX, mpz_t RY, mpz_t RZ,
                   vec_curve *curve,
                   mpz_t *X, mpz_t *Y, mpz_t *Z,
                   mpz_t *scalars,
                   size_t len)
{
  felem rx;
  felem ry;
  felem rz;

  felem *x = mpz_t_s_to_felems(X, len);
  felem *y = mpz_t_s_to_felems(Y, len);
  felem *z = mpz_t_s_to_felems(Z, len);

  vec_jsmul_nistp224_inner(rx, ry, rz,
                           curve,
                           x, y, z,
                           scalars,
                           len);
  felem_to_mpz_t(RX, rx);
  felem_to_mpz_t(RY, ry);
  felem_to_mpz_t(RZ, rz);

  free(x);
  free(y);
  free(z);
}

vec_jfmul_tab_ptr
vec_jfmul_precomp_nistp224(vec_curve *curve,
                           mpz_t X, mpz_t Y, mpz_t Z,
                           size_t len)
{
  felem x;
  felem y;
  felem z;
  vec_jfmul_tab_ptr ptr;

  ptr.nistp224 =
    (vec_jfmul_tab_nistp224_inner*)
    malloc(sizeof(vec_jfmul_tab_nistp224_inner));

  mpz_t_to_felem(x, X);
  mpz_t_to_felem(y, Y);
  mpz_t_to_felem(z, Z);

  vec_jfmul_init_nistp224_inner(ptr.nistp224, curve, len);

  vec_jfmul_prcmp_nistp224_inner(curve, ptr.nistp224, x, y, z);

  return ptr;
}

void
vec_jfmul_nistp224(mpz_t RX, mpz_t RY, mpz_t RZ,
                   vec_curve *curve,
                   vec_jfmul_tab_ptr ptr,
                   mpz_t scalar)
{
  felem rx;
  felem ry;
  felem rz;

  vec_jfmul_cmp_nistp224_inner(rx, ry, rz,
                               curve, ptr.nistp224,
                               scalar);

  felem_to_mpz_t(RX, rx);
  felem_to_mpz_t(RY, ry);
  felem_to_mpz_t(RZ, rz);
}

void
vec_jfmul_free_nistp224(vec_jfmul_tab_ptr ptr)
{
  vec_jfmul_clear_free_nistp224_inner(ptr.nistp224);
}

void
vec_jdbl_nistp224(vec_scratch_mpz_t scratch,
                  mpz_t X3, mpz_t Y3, mpz_t Z3,
                  vec_curve *curve,
                  mpz_t X1, mpz_t Y1, mpz_t Z1)
{
  felem x_in;
  felem y_in;
  felem z_in;
  felem x_out;
  felem y_out;
  felem z_out;

  VEC_UNUSED(curve);
  VEC_UNUSED(scratch);

  mpz_t_to_felem(x_in, X1);
  mpz_t_to_felem(y_in, Y1);
  mpz_t_to_felem(z_in, Z1);

  point_double(x_out, y_out, z_out, x_in, y_in, z_in);

  felem_to_mpz_t(X3, x_out);
  felem_to_mpz_t(Y3, y_out);
  felem_to_mpz_t(Z3, z_out);
}

void
vec_jadd_nistp224(vec_scratch_mpz_t scratch,
                  mpz_t X3, mpz_t Y3, mpz_t Z3,
                  vec_curve *curve,
                  mpz_t X1, mpz_t Y1, mpz_t Z1,
                  mpz_t X2, mpz_t Y2, mpz_t Z2)
{
  felem x1;
  felem y1;
  felem z1;
  felem x2;
  felem y2;
  felem z2;
  felem x3;
  felem y3;
  felem z3;

  VEC_UNUSED(scratch);
  VEC_UNUSED(curve);

  mpz_t_to_felem(x1, X1);
  mpz_t_to_felem(y1, Y1);
  mpz_t_to_felem(z1, Z1);
  mpz_t_to_felem(x2, X2);
  mpz_t_to_felem(y2, Y2);
  mpz_t_to_felem(z2, Z2);

  point_add(x3, y3, z3, x1, y1, z1, x2, y2, z2);

  felem_to_mpz_t(X3, x3);
  felem_to_mpz_t(Y3, y3);
  felem_to_mpz_t(Z3, z3);
}

/* These are timing routines and not tested beyond using them. */
/* LCOV_EXCL_START */

long
time_jdbl_nistp224(int test_time, mpz_t X, mpz_t Y)
{
  long i;
  int t;
  felem x;
  felem y;
  felem z;

  mpz_t Z;

  mpz_t_to_felem(x, X);
  mpz_t_to_felem(y, Y);

  mpz_init(Z);
  mpz_set_ui(Z, 1);
  mpz_t_to_felem(z, Z);
  mpz_clear(Z);

  t = clock();

  i = 0;
  do
    {
      point_double(x, y, z, x, y, z);
      i++;
    }
  while (!vec_done(t, test_time));

  return i;
}

long
time_jadd_nistp224(int test_time, mpz_t X, mpz_t Y)
{
  long i;
  int t;

  felem rx;
  felem ry;
  felem rz;

  felem x;
  felem y;
  felem z;

  mpz_t Z;

  mpz_t_to_felem(x, X);
  mpz_t_to_felem(y, Y);

  mpz_init(Z);
  mpz_set_ui(Z, 1);
  mpz_t_to_felem(z, Z);
  mpz_clear(Z);

  point_double(rx, ry, rz, x, y, z);

  t = clock();

  i = 0;
  do
    {
      point_add(rx, ry, rz, rx, ry, rz, x, y, z);
      i++;
    }
  while (!vec_done(t, test_time));

  return i;
}
/* LCOV_EXCL_STOP */
