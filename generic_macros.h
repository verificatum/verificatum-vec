
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

#include "undefine_macros.h"

#define SCRATCH(scratch) vec_scratch_mpz_t scratch

#define POSTFIX _generic_inner
#define TAB_POSTFIX _generic_inner

#define FIELD_ELEMENT mpz_t
#define FIELD_ELEMENT_INIT(x) mpz_init(x)
#define FIELD_ELEMENT_CLEAR(x) mpz_clear(x)
#define FIELD_ELEMENT_UNIT(x, y, z) \
  mpz_set_si(x, 0);        \
  mpz_set_si(y, 1);        \
  mpz_set_si(z, 0)
#define FIELD_ELEMENT_SET(rx, ry, rz, x, y, z)    \
  mpz_set(rx, x);                                 \
  mpz_set(ry, y);                                 \
  mpz_set(rz, z)

#define FIELD_ELEMENT_VAR mpz_t
#define FIELD_ELEMENT_VAR_INIT(x) mpz_init(x)
#define FIELD_ELEMENT_VAR_CLEAR(x) mpz_clear(x)
#define FIELD_ELEMENT_VAR_UNIT(x, y, z) \
  mpz_set_si(x, 0);                     \
  mpz_set_si(y, 1);                     \
  mpz_set_si(z, 0)
#define FIELD_ELEMENT_VAR_SET(rx, ry, rz, x, y, z) \
  mpz_set(rx, x);                                  \
  mpz_set(ry, y);                                  \
  mpz_set(rz, z)

#define FIELD_ELEMENT_CONTRACT(rx, ry, rz, x, y, z) \
  mpz_set(rx, x);                                   \
  mpz_set(ry, y);                                   \
  mpz_set(rz, z)                                    \

#define SCRATCH_INIT(scratch) vec_scratch_init_mpz_t(scratch)
#define SCRATCH_CLEAR(scratch) vec_scratch_clear_mpz_t(scratch)

#define ARRAY_MALLOC_INIT(len) vec_array_alloc_init(len)

#define ARRAY_CLEAR_FREE(array, len) vec_array_clear_free(array, len)

#define JDBL(scratch, rx, ry, rz, curve, x, y, z) \
  curve->jdbl(scratch,                            \
              rx, ry, rz,                         \
              curve,                              \
              x, y, z)

#define JDBL_VAR(scratch, rx, ry, rz, curve, x, y, z)  \
  curve->jdbl(scratch,                                 \
              rx, ry, rz,                              \
              curve,                                   \
              x, y, z)

#define JADD(scratch, rx, ry, rz, curve, x1, y1, z1, x2, y2, z2) \
  curve->jadd(scratch,                                           \
              rx, ry, rz,                                        \
              curve,                                             \
              x1, y1, z1,                                        \
              x2, y2, z2)

#define JADD_VAR(scratch, rx, ry, rz, curve, x1, y1, z1, x2, y2, z2)  \
  curve->jadd(scratch,                                                \
              rx, ry, rz,                                             \
              curve,                                                  \
              x1, y1, z1,                                             \
              x2, y2, z2)

#define CURVE vec_curve
