
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

#include <stdio.h>
#include <gmp.h>
#include "vec.h"

void
vec_jsmul_aff(mpz_t ropx, mpz_t ropy,
              vec_curve *curve,
              mpz_t *basesx, mpz_t *basesy,
              mpz_t *exponents,
              size_t len)
{
  size_t i;

  mpz_t ropz;
  mpz_t *basesz = vec_array_alloc_init(len);

  mpz_init(ropz);

  for (i = 0; i < len; i++) {
    vec_affj(basesx[i], basesy[i], basesz[i]);
  }

  curve->jsmul(ropx, ropy, ropz,
               curve,
               basesx, basesy, basesz,
               exponents,
               len);

  vec_jaff(ropx, ropy, ropz, curve);

  mpz_clear(ropz);

  vec_array_clear_free(basesz, len);
}
