
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
vec_smul_block_batch(mpz_t ropx, mpz_t ropy,
                     vec_curve *curve,
                     mpz_t *basesx, mpz_t *basesy,
                     mpz_t *scalars,
                     size_t len,
                     size_t block_width, size_t batch_len,
                     size_t max_scalar_bitlen)
{
  size_t i;
  vec_smul_tab table;
  mpz_t tmpx;
  mpz_t tmpy;

  vec_scratch_mpz_t scratch;

  mpz_init(tmpx);
  mpz_init(tmpy);

  vec_scratch_init_mpz_t(scratch);

  if (len < batch_len) {
    batch_len = len;
  }

  vec_smul_init(table, curve, batch_len, block_width);

  /* Initialize result to unit element. */
  mpz_set_si(ropx, -1);
  mpz_set_si(ropy, -1);

  for (i = 0; i < len; i += batch_len)
    {

      /* Last batch may be slightly shorter. */
      if (len - i < batch_len)
        {
          batch_len = len - i;

          vec_smul_clear(table);
          vec_smul_init(table, curve, batch_len, block_width);
        }

      /* Perform computation for batch */
      vec_smul_precomp(table, basesx, basesy);

      /* Compute batch. */
      vec_smul_table(tmpx, tmpy, table, scalars, max_scalar_bitlen);

      /* Add with result so far. */
      vec_add(scratch,
              ropx, ropy,
              curve,
              ropx, ropy,
              tmpx, tmpy);

      /* Move on to next batch. */
      basesx += batch_len;
      basesy += batch_len;
      scalars += batch_len;
    }

  vec_scratch_clear_mpz_t(scratch);

  mpz_clear(tmpy);
  mpz_clear(tmpx);

  vec_smul_clear(table);
}
