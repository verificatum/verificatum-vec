
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

void
vec_jaff(mpz_t X, mpz_t Y, mpz_t Z, vec_curve *curve) {

  mpz_t ZZ;

  if (mpz_cmp_si(Z, 0) == 0)
    {
      mpz_set_si(X, -1);
      mpz_set_si(Y, -1);
    }
  else
    {

      mpz_init(ZZ);

      mpz_invert(Z, Z, curve->modulus);

      mpz_mul(ZZ, Z, Z);
      mpz_mod(ZZ, ZZ, curve->modulus);

      mpz_mul(X, X, ZZ);
      mpz_mod(X, X, curve->modulus);

      mpz_mul(Z, ZZ, Z);
      mpz_mod(Z, Z, curve->modulus);

      mpz_mul(Y, Y, Z);
      mpz_mod(Y, Y, curve->modulus);

      mpz_clear(ZZ);
    }
}
