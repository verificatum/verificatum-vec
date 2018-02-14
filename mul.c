
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
vec_mul(mpz_t rx, mpz_t ry,
        vec_curve *curve,
        mpz_t x, mpz_t y,
        mpz_t scalar)
{

  int i;
  int bit_length;

  vec_scratch_mpz_t scratch;

  vec_scratch_init_mpz_t(scratch);

  /* Initialize with the unit element. */
  mpz_set_si(rx, -1);
  mpz_set_si(ry, -1);

  /* Determine bit length. */
  bit_length = (int)mpz_sizeinbase(scalar, 2);

  for (i = bit_length; i >= 0; i--)
    {

      /* Double. */
      vec_dbl(scratch,
              rx, ry,
              curve,
              rx, ry);

      /* Add. */
      if (mpz_tstbit(scalar, i))
        {

          vec_add(scratch,
                  rx, ry,
                  curve,
                  rx, ry,
                  x, y);
        }
    }

  vec_scratch_clear_mpz_t(scratch);
}
