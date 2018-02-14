
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
vec_jdbl_aff(vec_scratch_mpz_t scratch,
             mpz_t rx, mpz_t ry,
             vec_curve *curve,
             mpz_t x, mpz_t y) {
  mpz_t X1;
  mpz_t Y1;
  mpz_t Z1;
  mpz_t Z3;

  mpz_init(X1);
  mpz_init(Y1);
  mpz_init(Z1);
  mpz_init(Z3);

  mpz_set(X1, x);
  mpz_set(Y1, y);

  vec_affj(X1, Y1, Z1);

  curve->jdbl(scratch,
              rx, ry, Z3,
              curve,
              X1, Y1, Z1);

  vec_jaff(rx, ry, Z3, curve);

  mpz_clear(Z3);
  mpz_clear(Z1);
  mpz_clear(Y1);
  mpz_clear(X1);
}
