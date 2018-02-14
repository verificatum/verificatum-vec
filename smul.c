
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
vec_smul(mpz_t ropx, mpz_t ropy,
         vec_curve *curve,
         mpz_t *basesx, mpz_t *basesy,
         mpz_t *scalars,
         size_t len)
{
  size_t i;
  size_t bitlen;
  size_t max_scalar_bitlen;
  size_t batch_len = 100;      /* This is somewhat arbitrary, but it
                                  makes the amortized cost for squaring
                                  very small in comparison to the cost
                                  for multiplications. */
  size_t block_width;

  /* Compute the maximal bit length among the scalars. */
  max_scalar_bitlen = 0;
  for (i = 0; i < len; i++)
    {
      bitlen = mpz_sizeinbase(scalars[i], 2);
      if (bitlen > max_scalar_bitlen)
        {
          max_scalar_bitlen = bitlen;
        }
    }

  /* Determine a good block width. */
  block_width = vec_smul_block_width(max_scalar_bitlen, batch_len);

  vec_smul_block_batch(ropx, ropy,
                       curve,
                       basesx, basesy,
                       scalars,
                       len,
                       block_width, batch_len,
                       max_scalar_bitlen);
}
