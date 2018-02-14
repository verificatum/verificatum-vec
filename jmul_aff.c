
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
vec_jmul_aff(mpz_t rx, mpz_t ry,
             vec_curve *curve,
             mpz_t x, mpz_t y,
             mpz_t scalar) {

  mpz_t RZ;

  mpz_t X;
  mpz_t Y;
  mpz_t Z;

  mpz_init(RZ);

  mpz_init(X);
  mpz_init(Y);
  mpz_init(Z);

  mpz_set(X, x);
  mpz_set(Y, y);

  vec_affj(X, Y, Z);

  curve->jmul(rx, ry, RZ,
              curve,
              X, Y, Z,
              scalar);

  vec_jaff(rx, ry, RZ, curve);

  mpz_clear(Z);
  mpz_clear(Y);
  mpz_clear(X);

  mpz_clear(RZ);

}
