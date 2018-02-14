
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
#include <stdlib.h>
#include <gmp.h>

#include "vec.h"

vec_curve*
vec_curve_alloc()
{
  vec_curve *curve = (vec_curve *)malloc(sizeof(vec_curve));

  mpz_init(curve->modulus);
  mpz_init(curve->a);
  mpz_init(curve->b);
  mpz_init(curve->gx);
  mpz_init(curve->gy);
  mpz_init(curve->n);

  return curve;
}
