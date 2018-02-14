
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

#include <gmp.h>

#include "vec.h"

vec_curve *
vec_curve_get(mpz_t modulus, mpz_t a, mpz_t b,
                mpz_t gx, mpz_t gy, mpz_t n)
{
  vec_curve *curve = vec_curve_alloc();

  curve->name = NULL;

  mpz_set(curve->modulus, modulus);
  mpz_set(curve->a, a);
  mpz_set(curve->b, b);

  mpz_set(curve->gx, gx);
  mpz_set(curve->gy, gy);
  mpz_set(curve->n, n);

  curve->jdbl = vec_jdbl;
  curve->jadd = vec_jadd;
  curve->jmul = vec_jmul;
}
