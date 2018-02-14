
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
vec_smul_precomp(vec_smul_tab table, mpz_t *basesx, mpz_t *basesy)
{
  size_t i, j;        /* Index variables. */
  size_t block_width; /* Width of current subtable. */
  mpz_t *tx;          /* Temporary variable for subtable. */
  mpz_t *ty;          /* Temporary variable for subtable. */

  vec_scratch_mpz_t scratch;

  int mask;           /* Mask used for dynamic programming */
  int one_mask;       /* Mask containing a single non-zero bit. */

  vec_scratch_init_mpz_t(scratch);

  block_width = table->block_width;

  for (i = 0; i < table->tabs_len; i++)
    {
      /* Last block may have smaller width, but the width is never
         zero. */
      if (i == table->tabs_len - 1)
        {
          block_width = table->len - (table->tabs_len - 1) * block_width;
        }

      /* Current subtable. */
      tx = table->tabsx[i];
      ty = table->tabsy[i];

      /* Initialize current subtable with all trivial products. */
      mpz_set_si(tx[0], -1);
      mpz_set_si(ty[0], -1);

      mask = 1;
      for (j = 0; j < block_width; j++)
        {
          mpz_set(tx[mask], basesx[j]);
          mpz_set(ty[mask], basesy[j]);
          mask <<= 1;
        }

      /* Initialize current subtable with all non-trivial products. */
      for (mask = 1; mask < (1 << block_width); mask++)
        {
          one_mask = mask & (-mask);

          vec_add(scratch,
                  tx[mask], ty[mask],
                  table->curve,
                  tx[mask ^ one_mask], ty[mask ^ one_mask],
                  tx[one_mask], ty[one_mask]);
        }

      basesx += block_width;
      basesy += block_width;
    }

  vec_scratch_clear_mpz_t(scratch);
}
