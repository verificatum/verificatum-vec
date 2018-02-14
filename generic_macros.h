
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
