
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

#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gmp.h>

#include "vec.h"

#define DEFAULT_TEST_TIME 500
#define DEFAULT_SPEED_TIME 2000


/*
 * We use compiler flags that enforce that unused variables are
 * flagged as errors. Sometimes we need to declare an exception
 * explicitly state explicitly.
 */
#define VEC_UNUSED(x) ((void)(x))

/* LCOV_EXCL_START */
void
fail(char *str)
{
  fprintf(stderr, "FAILED!\n");
  fprintf(stderr, "%s", str);
  exit(-1);
}
/* LCOV_EXCL_STOP */

void
test_dbl_add(vec_curve *curve)
{

  int t;

  mpz_t tx;
  mpz_t ty;
  mpz_t rx1;
  mpz_t ry1;
  mpz_t rx2;
  mpz_t ry2;

  mpz_t Ox;
  mpz_t Oy;

  vec_scratch_mpz_t scratch;

  mpz_init(tx);
  mpz_init(ty);
  mpz_init(rx1);
  mpz_init(ry1);
  mpz_init(rx2);
  mpz_init(ry2);

  mpz_init(Ox);
  mpz_init(Oy);

  vec_scratch_init_mpz_t(scratch);

  mpz_set_si(Ox, -1);
  mpz_set_si(Oy, -1);

  t = clock();

  /* Test doubling unit element. */
  vec_dbl(scratch,
          rx1, ry1,
          curve,
          Ox, Oy);
  assert(vec_eq(rx1, ry1, Ox, Oy));

  /* Test adding unit element. */
  vec_add(scratch,
          rx1, ry1,
          curve,
          curve->gx, curve->gy,
          Ox, Oy);
  assert(vec_eq(rx1, ry1, curve->gx, curve->gy));

  vec_add(scratch,
          rx1, ry1,
          curve,
          Ox, Oy,
          curve->gx, curve->gy);
  assert(vec_eq(rx1, ry1, curve->gx, curve->gy));

    /* Test adding point to its negation. */
  mpz_sub(rx2, curve->modulus, curve->gy);
  vec_add(scratch,
          rx1, ry1,
          curve,
          curve->gx, rx2,
          curve->gx, curve->gy);
  assert(vec_eq(rx1, ry1, Ox, Oy));

  /* Test processing of general elements. */
  mpz_set(tx, curve->gx);
  mpz_set(ty, curve->gy);

  do
    {

      /* 4 * (tx, ty) */
      vec_dbl(scratch,
              rx1, ry1,
              curve,
              tx, ty);
      vec_dbl(scratch,
              rx1, ry1,
              curve,
              rx1, ry1);

      /* 2 * (tx, ty) + (tx, ty) + (tx, ty) */
      vec_dbl(scratch,
              rx2, ry2,
              curve,
              tx, ty);
      vec_add(scratch,
              rx2, ry2,
              curve,
              rx2, ry2,
              tx, ty);
      vec_add(scratch,
              rx2, ry2,
              curve,
              rx2, ry2,
              tx, ty);

      assert(vec_eq(rx1, ry1, rx2, ry2));

      mpz_set(tx, rx1);
      mpz_set(ty, ry1);
    }
  while (!vec_done(t, DEFAULT_TEST_TIME));

  vec_scratch_clear_mpz_t(scratch);

  mpz_clear(Ox);
  mpz_clear(Oy);

  mpz_clear(ry2);
  mpz_clear(rx2);
  mpz_clear(ry1);
  mpz_clear(rx1);
  mpz_clear(ty);
  mpz_clear(tx);
}

void
test_mul(vec_curve *curve) {

  int t;

  mpz_t rx1;
  mpz_t ry1;
  mpz_t rx2;
  mpz_t ry2;
  mpz_t rx3;
  mpz_t ry3;
  mpz_t rx4;
  mpz_t ry4;

  vec_scratch_mpz_t scratch;

  mpz_t Ox;
  mpz_t Oy;

  mpz_t scalar1;
  mpz_t scalar2;
  mpz_t scalar3;

  mpz_init(rx1);
  mpz_init(ry1);
  mpz_init(rx2);
  mpz_init(ry2);
  mpz_init(rx3);
  mpz_init(ry3);
  mpz_init(rx4);
  mpz_init(ry4);

  mpz_init(Ox);
  mpz_init(Oy);

  vec_scratch_init_mpz_t(scratch);

  mpz_init(scalar1);
  mpz_init(scalar2);
  mpz_init(scalar3);

  mpz_set_si(Ox, -1);
  mpz_set_si(Oy, -1);

  /* Generate "random" scalar1. */
  mpz_set_ui(scalar1, 1);
  mpz_mul_2exp(scalar1, scalar1, 1000);
  mpz_mod(scalar1, scalar1, curve->n);

  /* Test multiplying unit element. */
  vec_mul(rx1, ry1,
          curve,
          Ox, Oy,
          scalar1);
  assert(vec_eq(rx1, ry1, Ox, Oy));

  /* Test general multiplication. */
  mpz_set_ui(scalar2, 1);
  mpz_mul_2exp(scalar2, scalar2, 1000);
  mpz_mod(scalar2, scalar2, curve->n);

  t = clock();

  do
    {

      vec_mul(rx1, ry1,
              curve,
              curve->gx, curve->gy,
              scalar1);

      vec_mul(rx2, ry2,
              curve,
              curve->gx, curve->gy,
              scalar2);

      mpz_add(scalar3, scalar1, scalar2);
      vec_mul(rx3, ry3,
              curve,
              curve->gx, curve->gy,
              scalar3);

      vec_add(scratch,
              rx4, ry4,
              curve,
              rx1, ry1,
              rx2, ry2);

      assert(vec_eq(rx4, ry4, rx3, ry3));

      mpz_add(scalar1, scalar1, scalar2);
      mpz_mod(scalar1, scalar1, curve->n);

      mpz_mul(scalar2, scalar2, scalar2);
      mpz_mod(scalar2, scalar2, curve->n);

    }
  while (!vec_done(t, DEFAULT_TEST_TIME));

  mpz_clear(scalar2);
  mpz_clear(scalar1);

  vec_scratch_clear_mpz_t(scratch);

  mpz_clear(Ox);
  mpz_clear(Oy);

  mpz_clear(ry4);
  mpz_clear(rx4);
  mpz_clear(ry3);
  mpz_clear(rx3);
  mpz_clear(ry2);
  mpz_clear(rx2);
  mpz_clear(ry1);
  mpz_clear(rx1);
}

void
test_smul(vec_curve *curve) {

  int t;
  size_t len;
  size_t i;

  mpz_t rx1;
  mpz_t ry1;

  mpz_t rx2;
  mpz_t ry2;

  mpz_t *basesx;
  mpz_t *basesy;
  mpz_t *scalars;

  mpz_t scalar;

  vec_scratch_mpz_t scratch;

  mpz_init(rx1);
  mpz_init(ry1);

  mpz_init(rx2);
  mpz_init(ry2);

  mpz_init(scalar);

  vec_scratch_init_mpz_t(scratch);

  mpz_set_ui(scalar, 1);
  mpz_mul_2exp(scalar, scalar, 100000);
  mpz_mod(scalar, scalar, curve->n);

  len = 1;

  t = clock();

  do
    {

      /* Generate "random" bases and scalars. */
      basesx = vec_array_alloc_init(len);
      basesy = vec_array_alloc_init(len);
      scalars = vec_array_alloc_init(len);

      for (i = 0; i < len; i++) {

        vec_mul(basesx[i], basesy[i],
                curve,
                curve->gx, curve->gy,
                scalar);

        mpz_mul(scalar, scalar, scalar);
        mpz_mod(scalar, scalar, curve->n);
        mpz_set(scalars[i], scalar);
      }

      /* Compute simultaneous multiplication naively. */
      mpz_set_si(rx1, -1);
      mpz_set_si(ry1, -1);

      for (i = 0; i < len; i++) {

        vec_mul(rx2, ry2,
                curve,
                basesx[i], basesy[i],
                scalars[i]);

        vec_add(scratch,
                rx1, ry1,
                curve,
                rx1, ry1,
                rx2, ry2);
      }

      vec_smul(rx2, ry2,
               curve,
               basesx, basesy,
               scalars,
               len);

      assert(vec_eq(rx2, ry2, rx1, ry1));


      vec_array_clear_free(scalars, len);
      vec_array_clear_free(basesy, len);
      vec_array_clear_free(basesx, len);

      len <<= 1;
    }
  while (!vec_done(t, DEFAULT_TEST_TIME));

  vec_scratch_clear_mpz_t(scratch);

  mpz_clear(scalar);

  mpz_clear(ry2);
  mpz_clear(rx2);
  mpz_clear(ry1);
  mpz_clear(rx1);

}

void
test_jdbl(vec_curve *curve)
{

  int t;

  mpz_t rx1;
  mpz_t ry1;
  mpz_t rx2;
  mpz_t ry2;

  vec_scratch_mpz_t scratch;

  mpz_t Ox;
  mpz_t Oy;

  mpz_init(rx1);
  mpz_init(ry1);
  mpz_init(rx2);
  mpz_init(ry2);

  mpz_init(Ox);
  mpz_init(Oy);

  vec_scratch_init_mpz_t(scratch);

  t = clock();

  mpz_set_si(Ox, -1);
  mpz_set_si(Oy, -1);

  /* Test processing of unit element. */
  vec_jdbl_aff(scratch,
               rx1, ry1,
               curve,
               Ox, Oy);
  assert(vec_eq(rx1, ry1, Ox, Oy));

  /* Test processing of general elements. */
  mpz_set(rx1, curve->gx);
  mpz_set(ry1, curve->gy);

  mpz_set(rx2, curve->gx);
  mpz_set(ry2, curve->gy);

  do
    {

      vec_dbl(scratch,
              rx1, ry1,
              curve,
              rx1, ry1);

      vec_jdbl_aff(scratch,
                   rx2, ry2,
                   curve,
                   rx2, ry2);

      assert(vec_eq(rx1, ry1, rx2, ry2));
    }
  while (!vec_done(t, DEFAULT_TEST_TIME));

  vec_scratch_init_mpz_t(scratch);

  mpz_clear(Oy);
  mpz_clear(Ox);

  mpz_clear(ry2);
  mpz_clear(rx2);
  mpz_clear(ry1);
  mpz_clear(rx1);
}

void
test_jadd(vec_curve *curve)
{

  int t;

  mpz_t rx1;
  mpz_t ry1;
  mpz_t rx2;
  mpz_t ry2;

  vec_scratch_mpz_t scratch;

  mpz_t Ox;
  mpz_t Oy;

  mpz_t bx;
  mpz_t by;

  mpz_init(rx1);
  mpz_init(ry1);
  mpz_init(rx2);
  mpz_init(ry2);

  mpz_init(Ox);
  mpz_init(Oy);

  mpz_init(bx);
  mpz_init(by);

  vec_scratch_init_mpz_t(scratch);

  t = clock();

  mpz_set_si(Ox, -1);
  mpz_set_si(Oy, -1);

  /* Test processing of unit element. */
  vec_jadd_aff(scratch,
               rx1, ry1,
               curve,
               Ox, Oy,
               Ox, Oy);
  assert(vec_eq(rx1, ry1, Ox, Oy));

  vec_jadd_aff(scratch,
               rx1, ry1,
               curve,
               curve->gx, curve->gy,
               Ox, Oy);
  assert(vec_eq(rx1, ry1, curve->gx, curve->gy));

  vec_jadd_aff(scratch,
               rx1, ry1,
               curve,
               Ox, Oy,
               curve->gx, curve->gy);
  assert(vec_eq(rx1, ry1, curve->gx, curve->gy));

  /* Test adding point to its negation. */
  mpz_sub(rx2, curve->modulus, curve->gy);
  vec_jadd_aff(scratch,
               rx1, ry1,
               curve,
               curve->gx, rx2,
               curve->gx, curve->gy);
  assert(vec_eq(rx1, ry1, Ox, Oy));

  /* Test processing of general elements. */
  mpz_set(rx1, curve->gx);
  mpz_set(ry1, curve->gy);

  mpz_set(rx2, curve->gx);
  mpz_set(ry2, curve->gy);

  mpz_set(bx, curve->gx);
  mpz_set(by, curve->gy);

  do
    {

      vec_add(scratch,
              rx1, ry1,
              curve,
              rx1, ry1,
              bx, by);

      vec_jadd_aff(scratch,
                   rx2, ry2,
                   curve,
                   rx2, ry2,
                   bx, by);

      vec_dbl(scratch,
              bx, by,
              curve,
              bx, by);
      vec_add(scratch,
              bx, by,
              curve,
              bx, by,
              curve->gx, curve->gy);

      assert(vec_eq(rx1, ry1, rx2, ry2));
    }
  while (!vec_done(t, DEFAULT_TEST_TIME));

  vec_scratch_clear_mpz_t(scratch);

  mpz_clear(by);
  mpz_clear(bx);

  mpz_clear(Oy);
  mpz_clear(Ox);

  mpz_clear(ry2);
  mpz_clear(rx2);
  mpz_clear(ry1);
  mpz_clear(rx1);
}

void
test_jmul(vec_curve *curve)
{
  int t;

  mpz_t rx1;
  mpz_t ry1;
  mpz_t rx2;
  mpz_t ry2;

  vec_scratch_mpz_t scratch;

  mpz_t Ox;
  mpz_t Oy;

  mpz_t bx;
  mpz_t by;

  mpz_t scalar;

  mpz_init(rx1);
  mpz_init(ry1);
  mpz_init(rx2);
  mpz_init(ry2);

  vec_scratch_init_mpz_t(scratch);

  mpz_init(Ox);
  mpz_init(Oy);

  mpz_init(bx);
  mpz_init(by);

  mpz_init(scalar);

  t = clock();

  mpz_set_si(Ox, -1);
  mpz_set_si(Oy, -1);

  /* Test multiplication with unit element. */
  mpz_set_ui(scalar, 1);
  mpz_mul_2exp(scalar, scalar, 100000);
  mpz_mod(scalar, scalar, curve->n);

  vec_jmul_aff(rx1, ry1,
               curve,
               Ox, Oy,
               scalar);
  assert(vec_eq(rx1, ry1, Ox, Oy));

  mpz_set(bx, curve->gx);
  mpz_set(by, curve->gy);

  /* Test general multiplication. */
  do
    {

      vec_mul(rx1, ry1,
              curve,
              bx, by,
              scalar);

      vec_jmul_aff(rx2, ry2,
                   curve,
                   bx, by,
                   scalar);

      assert(vec_eq(rx1, ry1, rx2, ry2));

      mpz_set(bx, rx1);
      mpz_set(by, ry1);

      mpz_mul(scalar, scalar, scalar);
      mpz_mod(scalar, scalar, curve->n);

    }
  while (!vec_done(t, DEFAULT_TEST_TIME));

  mpz_clear(bx);
  mpz_clear(by);

  mpz_clear(Ox);
  mpz_clear(Oy);

  vec_scratch_clear_mpz_t(scratch);

  mpz_clear(ry2);
  mpz_clear(rx2);
  mpz_clear(ry1);
  mpz_clear(rx1);
}

void
test_jsmul(vec_curve *curve)
{

  int t;
  size_t len;
  size_t i;

  mpz_t rx1;
  mpz_t ry1;

  mpz_t rx2;
  mpz_t ry2;

  mpz_t *basesx;
  mpz_t *basesy;
  mpz_t *scalars;

  mpz_t scalar;

  mpz_init(rx1);
  mpz_init(ry1);

  mpz_init(rx2);
  mpz_init(ry2);

  mpz_init(scalar);

  mpz_set_ui(scalar, 1);
  mpz_mul_2exp(scalar, scalar, 100000);
  mpz_mod(scalar, scalar, curve->n);

  len = 1;

  t = clock();
  do
    {

      /* Generate "random" bases and scalars. */
      basesx = vec_array_alloc_init(len);
      basesy = vec_array_alloc_init(len);
      scalars = vec_array_alloc_init(len);

      for (i = 0; i < len; i++) {

        vec_mul(basesx[i], basesy[i],
                curve,
                curve->gx, curve->gy,
                scalar);

        mpz_mul(scalar, scalar, scalar);
        mpz_mod(scalar, scalar, curve->n);
        mpz_set(scalars[i], scalar);
      }

      /* Compute simultaneous multiplication affinely. */
      vec_smul(rx1, ry1,
               curve,
               basesx, basesy,
               scalars,
               len);

      /* Compute simultaneous multiplication Jacobi. */
      vec_jsmul_aff(rx2, ry2,
                    curve,
                    basesx, basesy,
                    scalars,
                    len);

      assert(vec_eq(rx2, ry2, rx1, ry1));

      vec_array_clear_free(scalars, len);
      vec_array_clear_free(basesy, len);
      vec_array_clear_free(basesx, len);

      len <<= 1;
    }
  while (!vec_done(t, DEFAULT_TEST_TIME));

  mpz_clear(scalar);

  mpz_clear(ry2);
  mpz_clear(rx2);
  mpz_clear(ry1);
  mpz_clear(rx1);

}

void
test_jfmul(vec_curve *curve)
{
  int t;

  mpz_t rx1;
  mpz_t ry1;
  mpz_t rx2;
  mpz_t ry2;

  vec_scratch_mpz_t scratch;
  vec_jfmul_tab_ptr table_ptr;

  mpz_t scalar;

  mpz_init(rx1);
  mpz_init(ry1);
  mpz_init(rx2);
  mpz_init(ry2);

  vec_scratch_init_mpz_t(scratch);

  mpz_init(scalar);

  t = clock();

  table_ptr = vec_jfmul_precomp_aff(curve, curve->gx, curve->gy, 1000);

  /* Test multiplication with unit element. */
  mpz_set_ui(scalar, 0);

  vec_mul(rx1, ry1,
          curve,
          curve->gx, curve->gy,
          scalar);

  vec_jfmul_aff(rx2, ry2,
                curve,
                table_ptr,
                scalar);

  assert(vec_eq(rx1, ry1, rx2, ry2));

  /* Test general multiplication. */

  mpz_set_ui(scalar, 1);
  mpz_mul_2exp(scalar, scalar, 100000);
  mpz_mod(scalar, scalar, curve->n);

  do
    {

      vec_mul(rx1, ry1,
              curve,
              curve->gx, curve->gy,
              scalar);

      vec_jfmul_aff(rx2, ry2,
                    curve,
                    table_ptr,
                    scalar);

      assert(vec_eq(rx1, ry1, rx2, ry2));

      mpz_mul(scalar, scalar, scalar);
      mpz_mod(scalar, scalar, curve->n);

    }
  while (!vec_done(t, DEFAULT_TEST_TIME));

  vec_jfmul_free_aff(curve, table_ptr);

  mpz_clear(scalar);

  vec_scratch_clear_mpz_t(scratch);

  mpz_clear(ry2);
  mpz_clear(rx2);
  mpz_clear(ry1);
  mpz_clear(rx1);
}

void
test_sqrt(mpz_t p) {

  int t;

  mpz_t a;
  mpz_t z;
  mpz_t res;

  mpz_init(a);
  mpz_init(z);
  mpz_init(res);

  t = clock();

  mpz_set_si(a, 1);
  mpz_mul_2exp(a, a, 1000);
  mpz_mod(a, a, p);

  mpz_set_si(z, 2);
  while (mpz_legendre(z, p) == 1)
    {
      mpz_add_ui(z, z, 1);
    }

  mpz_mul(a, a, a);
  mpz_mod(a, a, p);

  do
    {

      vec_sqrt(res, a, p);

      mpz_mul(res, res, res);
      mpz_mod(res, res, p);

      assert(mpz_cmp(res, a) == 0);

      /* Randomize a new square */
      mpz_powm(a, a, a, p);
      mpz_mul(a, a, a);
      mpz_mod(a, a, p);

    }
  while (!vec_done(t, DEFAULT_TEST_TIME));

  mpz_clear(res);
  mpz_clear(z);
  mpz_clear(a);
}


/* These are timing routines and not tested beyond using them. */
/* LCOV_EXCL_START */

void
print_time()
{
  time_t timer;
  char buffer[25];
  struct tm* tm_info;

  time(&timer);
  tm_info = localtime(&timer);

  strftime(buffer, 25, "%Y:%m:%d%H:%M:%S\n", tm_info);
  fprintf(stderr, "%s", buffer);
}

void
print_doublings(char *prefix, long ops)
{
  printf("%12ld %s doublings\n", ops, prefix);
  fflush(stdout);
}
void
print_additions(char *prefix, long ops)
{
  printf("%12ld %s additions\n", ops, prefix);
  fflush(stdout);
}
void
print_multiplications(char *prefix, long ops)
{
  printf("%12ld %s multiplications\n", ops, prefix);
  fflush(stdout);
}

void
print_test(char *str)
{
  printf("%s\n", str);
  fflush(stdout);
}

long
time_dbl(vec_curve *curve, jdbl_func jdbl, long millisecs)
{

  long i;
  int t;

  vec_scratch_mpz_t scratch;

  mpz_t x;
  mpz_t y;

  vec_scratch_init_mpz_t(scratch);

  mpz_init(x);
  mpz_init(y);

  mpz_set(x, curve->gx);
  mpz_set(y, curve->gy);

  t = clock();

  i = 0;
  do
    {
      if (jdbl != NULL)
        {
          vec_jdbl_aff(scratch,
                       x, y,
                       curve,
                       x, y);
        }
      else
        {
          vec_dbl(scratch,
                  x, y,
                  curve,
                  x, y);
        }
      i++;
    }
  while (!vec_done(t, millisecs));

  mpz_clear(y);
  mpz_clear(x);

  vec_scratch_clear_mpz_t(scratch);

  return i;
}

long
time_jdbl(vec_curve *curve, long millisecs)
{
  long i;
  int t;

  vec_scratch_mpz_t scratch;

  mpz_t GX;
  mpz_t GY;
  mpz_t GZ;

  mpz_t X;
  mpz_t Y;
  mpz_t Z;

  vec_scratch_init_mpz_t(scratch);

  mpz_init(GX);
  mpz_init(GY);
  mpz_init(GZ);

  mpz_init(X);
  mpz_init(Y);
  mpz_init(Z);

  mpz_set(GX, curve->gx);
  mpz_set(GY, curve->gy);
  mpz_set_si(GZ, 1);

  /* Make sure we measure with GZ != 1. */
  vec_jdbl_generic(scratch,
                   X, Y, Z,
                   curve,
                   GX, GY, GZ);

  vec_jadd_generic(scratch,
                   X, Y, Z,
                   curve,
                   X, Y, Z,
                   GX, GY, GZ);

  t = clock();

  i = 0;
  do
    {
      curve->jdbl(scratch,
                  X, Y, Z,
                  curve,
                  X, Y, Z);
      i++;
    }
  while (!vec_done(t, millisecs));

  mpz_clear(Z);
  mpz_clear(Y);
  mpz_clear(X);

  mpz_clear(GZ);
  mpz_clear(GY);
  mpz_clear(GX);

  vec_scratch_clear_mpz_t(scratch);

  return i;
}

long
time_add(vec_curve *curve, jadd_func jadd, long millisecs)
{

  long i;
  int t;

  vec_scratch_mpz_t scratch;

  mpz_t x;
  mpz_t y;

  vec_scratch_init_mpz_t(scratch);

  VEC_UNUSED(millisecs);

  mpz_init(x);
  mpz_init(y);

  mpz_set(x, curve->gx);
  mpz_set(y, curve->gy);

  t = clock();

  i = 0;
  do
    {
      if (jadd != NULL)
        {
          vec_jadd_aff(scratch,
                       x, y,
                       curve,
                       x, y,
                       curve->gx, curve->gy);
        }
      else
        {
          vec_add(scratch,
                  x, y,
                  curve,
                  x, y,
                  curve->gx, curve->gy);
        }
      i++;
    }
  while (!vec_done(t, millisecs));

  mpz_clear(y);
  mpz_clear(x);

  vec_scratch_clear_mpz_t(scratch);

  return i;
}

long
time_jadd(vec_curve *curve, long millisecs)
{
  long i;
  int t;

  vec_scratch_mpz_t scratch;

  mpz_t X;
  mpz_t Y;
  mpz_t Z;

  mpz_t GX;
  mpz_t GY;
  mpz_t GZ;

  vec_scratch_init_mpz_t(scratch);

  mpz_init(X);
  mpz_init(Y);
  mpz_init(Z);

  mpz_init(GX);
  mpz_init(GY);
  mpz_init(GZ);

  mpz_set(GX, curve->gx);
  mpz_set(GY, curve->gy);
  mpz_set_si(GZ, 1);

  /* Make sure we measure with GZ != 1. */
  vec_jdbl_generic(scratch,
                   GX, GY, GZ,
                   curve,
                   GX, GY, GZ);

  vec_jadd_generic(scratch,
                   GX, GY, GZ,
                   curve,
                   GX, GY, GZ,
                   GX, GY, GZ);

  mpz_set_si(X, 0);
  mpz_set_si(Y, 1);
  mpz_set_si(Z, 0);

  t = clock();

  i = 0;
  do
    {

      curve->jadd(scratch,
                  X, Y, Z,
                  curve,
                  X, Y, Z,
                  GX, GY, GZ);
      i++;
    }
  while (!vec_done(t, millisecs));

  mpz_clear(Z);
  mpz_clear(Y);
  mpz_clear(X);

  mpz_clear(GZ);
  mpz_clear(GY);
  mpz_clear(GX);

  vec_scratch_clear_mpz_t(scratch);

  return i;
}

long
time_mul(vec_curve *curve, jmul_func jmul, long millisecs)
{

  long i;
  int t;

  vec_scratch_mpz_t scratch;

  mpz_t x;
  mpz_t y;
  mpz_t bx;
  mpz_t by;

  mpz_t scalar;

  vec_scratch_init_mpz_t(scratch);

  mpz_init(x);
  mpz_init(y);
  mpz_init(bx);
  mpz_init(by);

  mpz_init(scalar);

  mpz_set(x, curve->gx);
  mpz_set(y, curve->gy);

  /* Test general multiplication. */
  mpz_set_ui(scalar, 1);
  mpz_mul_2exp(scalar, scalar, 123456);
  mpz_mod(scalar, scalar, curve->n);

  mpz_set(bx, curve->gx);
  mpz_set(by, curve->gy);

  t = clock();

  i = 0;
  do
    {
      if (jmul != NULL)
        {
          vec_jmul_aff(x, y,
                       curve,
                       bx, by,
                       scalar);
        }
      else
        {
          vec_mul(x, y,
                  curve,
                  bx, by,
                  scalar);
        }

      mpz_set(bx, x);
      mpz_set(by, y);

      mpz_mul_si(scalar, scalar, 2);
      mpz_mod(scalar, scalar, curve->n);

      i++;
    }
  while (!vec_done(t, millisecs));

  mpz_clear(scalar);

  mpz_clear(by);
  mpz_clear(bx);
  mpz_clear(y);
  mpz_clear(x);

  vec_scratch_clear_mpz_t(scratch);

  return i;
}

long
time_jfmul(vec_curve *curve, long millisecs)
{
  int t;
  long i;

  mpz_t rx;
  mpz_t ry;

  vec_scratch_mpz_t scratch;
  vec_jfmul_tab_ptr table_ptr;

  mpz_t scalar;

  VEC_UNUSED(millisecs);

  mpz_init(rx);
  mpz_init(ry);

  vec_scratch_init_mpz_t(scratch);

  mpz_init(scalar);

  /* Test general multiplication. */

  mpz_set_ui(scalar, 1);
  mpz_mul_2exp(scalar, scalar, 100000);
  mpz_mod(scalar, scalar, curve->n);

  t = clock();

  table_ptr = vec_jfmul_precomp_aff(curve, curve->gx, curve->gy, 400000);

  i = 0;
  do
    {

      vec_jfmul_aff(rx, ry,
                    curve,
                    table_ptr,
                    scalar);

      mpz_mul(scalar, scalar, scalar);
      mpz_mod(scalar, scalar, curve->n);

      i++;
    }
  while (!vec_done(t, DEFAULT_SPEED_TIME));

  curve->jfmul_free(table_ptr);

  mpz_clear(scalar);

  vec_scratch_clear_mpz_t(scratch);

  mpz_clear(ry);
  mpz_clear(rx);

  return i;
}

long
time_jmul(vec_curve *curve, long millisecs)
{
  long i;
  int t;

  mpz_t X;
  mpz_t Y;
  mpz_t Z;

  mpz_t BX;
  mpz_t BY;
  mpz_t BZ;

  mpz_t scalar;

  VEC_UNUSED(millisecs);

  mpz_init(X);
  mpz_init(Y);
  mpz_init(Z);

  mpz_init(BX);
  mpz_init(BY);
  mpz_init(BZ);

  mpz_init(scalar);

  /* Test general multiplication. */
  mpz_set_ui(scalar, 1);
  mpz_mul_2exp(scalar, scalar, 123456);
  mpz_mod(scalar, scalar, curve->n);

  mpz_set(BX, curve->gx);
  mpz_set(BY, curve->gy);
  mpz_set_si(BZ, 1);

  t = clock();

  i = 0;
  do
    {

      curve->jmul(X, Y, Z,
                  curve,
                  BX, BY, BZ,
                  scalar);

      mpz_set(BX, X);
      mpz_set(BY, Y);
      mpz_set(BZ, Z);

      mpz_mul_si(scalar, scalar, 2);
      mpz_mod(scalar, scalar, curve->n);

      i++;
    }
  while (!vec_done(t, millisecs));

  mpz_clear(scalar);

  mpz_clear(Z);
  mpz_clear(Y);
  mpz_clear(X);

  mpz_clear(BZ);
  mpz_clear(BY);
  mpz_clear(BX);

  return i;
}
/* LCOV_EXCL_STOP */

void
test_curve(char *name)
{

  vec_curve *curve = vec_curve_get_named(name, 0);

  if (curve == NULL)
    {
      /* LCOV_EXCL_START */
      fprintf(stderr, "Unknown curve name!\n");
      exit(1);
      /* LCOV_EXCL_STOP */
    }

  printf("\nTesting curve: %s (%ld ms/function)\n", curve->name,
         (long)DEFAULT_TEST_TIME);
  printf("----------------------------------------------------------------\n");

  print_test("Sqrt (solving quadratic equations)");
  test_sqrt(curve->modulus);

  print_test("Affine doubling and adding");
  test_dbl_add(curve);

  print_test("Affine multiplication");
  test_mul(curve);

  print_test("Affine simultaneous multiplication");
  test_smul(curve);

  print_test("Jacobi doubling");
  test_jdbl(curve);

  print_test("Jacobi adding");
  test_jadd(curve);

  print_test("Jacobi sliding-window multiplication");
  test_jmul(curve);

  print_test("Jacobi simultaneous multiplication");
  test_jsmul(curve);

  print_test("Jacobi fixed-basis multiplication");
  test_jfmul(curve);

  vec_curve_free(curve);

  curve = vec_curve_get_named(name, 1);

  if ((curve->jdbl != vec_jdbl_generic
       && curve->jdbl != vec_jdbl_a_eq_neg3_generic)
      || (curve->jadd != vec_jadd_generic
          && curve->jadd != vec_jadd_a_eq_neg3_generic)
      || (curve->jmul != vec_jmulsw_generic
          && curve->jmul != vec_jmulsw_a_eq_neg3_generic)
      || (curve->jsmul != vec_jsmul_generic
          && curve->jsmul != vec_jsmul_a_eq_neg3_generic))
    {
      printf("\nTesting optimized code for this curve.\n\n");
    }

  if (curve->jdbl != vec_jdbl_generic
      && curve->jdbl != vec_jdbl_a_eq_neg3_generic)
    {
      print_test("Jacobi doubling");
      test_jdbl(curve);
    }
  if (curve->jadd != vec_jadd_generic)
    {
      print_test("Jacobi adding");
      test_jadd(curve);
    }
  if (curve->jmul != vec_jmulsw_generic
      && curve->jmul != vec_jmulsw_a_eq_neg3_generic)
    {
      print_test("Jacobi sliding-window multiplication");
      test_jmul(curve);
    }
  if (curve->jsmul != vec_jsmul_generic)
    {
      print_test("Jacobi simultaneous multiplication");
      test_jsmul(curve);
    }
  if (curve->jfmul != vec_jfmul_generic)
    {
      print_test("Jacobi fixed-basis multiplication");
      test_jfmul(curve);
    }

  vec_curve_free(curve);
}

/* These are timing routines and not tested beyond using them. */
/* LCOV_EXCL_START */
void
time_curve(char *name, long millisecs)
{

  vec_curve *curve = vec_curve_get_named(name, 0);

  if (curve == NULL)
    {
      fprintf(stderr, "Unknown curve name!\n");
      exit(1);
    }

  printf("\nTiming curve: %s (%ld ms/function)\n", curve->name, millisecs);
  printf("----------------------------------------------------------------\n");

  print_doublings("Affine", time_dbl(curve, NULL, millisecs));
  print_additions("Affine", time_add(curve, NULL, millisecs));
  print_multiplications("Affine", time_mul(curve, NULL, millisecs));

  /* Jacobi. */
  print_doublings("Jacobi", time_jdbl(curve, millisecs));
  print_additions("Jacobi", time_jadd(curve, millisecs));
  print_multiplications("Jacobi sliding window", time_jmul(curve, millisecs));
  print_multiplications("Affined Jacobi sliding window",
                        time_mul(curve, curve->jmul, millisecs));
  print_multiplications("Affined Jacobi fixed-basis",
                        time_jfmul(curve, millisecs));

  vec_curve_free(curve);

  curve = vec_curve_get_named(name, 1);

  if (curve->jdbl != vec_jdbl_generic
      || curve->jadd != vec_jadd_generic
      || curve->jmul != vec_jmulsw_generic
      || curve->jsmul != vec_jsmul_generic)
    {

      printf("\nTiming optimized code for this curve.\n\n");

      if (curve->jdbl_timer != NULL)
        {
          print_doublings("Raw Jacobi",
                          curve->jdbl_timer(millisecs, curve->gx, curve->gy));
        }
      if (curve->jadd_timer != NULL)
        {
          print_additions("Raw Jacobi",
                          curve->jadd_timer(millisecs, curve->gx, curve->gy));
        }

      if (curve->jdbl != vec_jdbl_generic
          && curve->jdbl != vec_jdbl_a_eq_neg3_generic)
        {
          print_doublings("Jacobi", time_jdbl(curve, millisecs));
        }
      if (curve->jadd != vec_jadd_generic)
        {
          print_additions("Jacobi", time_jadd(curve, millisecs));
        }
      if (curve->jmul != vec_jmulsw_generic)
        {
          print_multiplications("Jacobi sliding window",
                                time_jmul(curve, millisecs));
          print_multiplications("Affined Jacobi sliding window",
                                time_mul(curve, curve->jmul, millisecs));
        }
      if (curve->jdbl != vec_jdbl_generic
          && curve->jdbl != vec_jdbl_a_eq_neg3_generic)
        {
          print_multiplications("Affined Jacobi fixed-basis",
                                time_jfmul(curve, millisecs));
        }
    }



  vec_curve_free(curve);
}

void
usage(char *command_name) {
  printf("Usage: %s check|speed [name ...]\n", command_name);
  exit(0);
}
/* LCOV_EXCL_STOP */

int
main(int args, char *argv[]) {

  int i;
  char *name = NULL;
  int test = 0;
  long millisecs = DEFAULT_SPEED_TIME;

  /* Command line */
  /* LCOV_EXCL_START */
  if (args == 1)
    {
      usage(argv[0]);
    }

  if (strcmp(argv[1], "check") == 0)
    {
      test = 1;
    }
  else if (strcmp(argv[1], "speed") == 0)
    {
      test = 0;
    }
  else
    {
      usage(argv[0]);
    }
  /* LCOV_EXCL_STOP */

  if (test == 0)
    {

      /* These are timing routines and not tested beyond using them. */
      /* LCOV_EXCL_START */

      printf("\n================================================================\n");
      printf("\n          BENCHMARKS FOR VEC\n\n");
      printf("You need to consult the code understand exactly what is\n");
      printf("measured before drawing any conclusions, but the benchmarks\n");
      printf("are fairly self explanatory.\n");
      printf("\n");
      printf("The default code is written by Douglas Wikstrom directly on\n");
      printf("top of the The GNU Multiple Precision Arithmetic Library.\n");
      printf("\n");
      printf("The optimized code is copied from the OpenSSL project\n");
      printf("http://www.openssl.org and written by Emilia Käsper and David\n");
      printf("Langley, inspired by code written by Dan Bernstein.\n");
      printf("The OpenSSL code is licensed under the Apache License 2.0.\n\n");
      printf("================================================================\n");
      /* LCOV_EXCL_START */

    } else {

    printf("\n================================================================\n");
    printf("\n      TESTS FOR VEC\n\n");
    printf("The default code is written by Douglas Wikstrom directly on\n");
    printf("of the The GNU Multiple Precision Arithmetic Library.\n");
    printf("\n");
    printf("The optimized code is copied from the OpenSSL project\n");
    printf("http://www.openssl.org and written by Emilia Käsper and David\n");
    printf("Langley, inspired by code written by Dan Bernstein.\n");
    printf("The OpenSSL code is licensed under the Apache License 2.0.\n\n");
    printf("================================================================\n");
  }

  if (args > 2)
    {
      /* Only used for manual command line use. */
      /* LCOV_EXCL_START */
      for (i = 2; i < args; i++)
        {
          name = argv[i];
          if (test)
            {
              test_curve(name);
            }
          else
            {
              time_curve(name, millisecs);
            }
        }
      /* LCOV_EXCL_START */
    }
  else
    {
      i = 0;

      while ((name = vec_curve_get_name(i)) != NULL)
        {
          if (test)
            {
              test_curve(name);
            }
          else
            {
              /* These are timing routines and not tested beyond using them. */
              /* LCOV_EXCL_START */
              time_curve(name, millisecs);
              /* LCOV_EXCL_START */
            }

          i++;
        }

      if (i != vec_curve_number_of_names()) {
        fail("Reported number of curve names is false!");
      }
    }

  return 0;
}
