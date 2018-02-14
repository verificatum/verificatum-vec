
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

void
vec_smul_init(vec_smul_tab table,
              vec_curve *curve,
              size_t len, size_t block_width)
{
  size_t i, j;     /* Index parameters. */
  size_t tab_len;  /* Size of a subtable. */
  mpz_t *tx;       /* Temporary variable for subtable. */
  mpz_t *ty;       /* Temporary variable for subtable. */

  /* Initialize the curve parameters of the table. */
  table->curve = curve;

  /* Initialize length parameters. */
  table->len = len;
  table->block_width = block_width;
  if (len < block_width) {
    table->block_width = len;
  }
  table->tabs_len = (len + block_width - 1) / block_width;

  /* Allocate and initialize space for pointers to tables. */
  table->tabsx = (mpz_t **)malloc(table->tabs_len * sizeof(mpz_t *));
  table->tabsy = (mpz_t **)malloc(table->tabs_len * sizeof(mpz_t *));

  tab_len = 1 << block_width;
  for (i = 0; i < table->tabs_len; i++)
    {

      /* Last block may be more narrow than the other, but it is never
         of size zero. */
      if (i == table->tabs_len - 1
          && len - (table->tabs_len - 1) * block_width < block_width)
        {
          block_width = len - (table->tabs_len - 1) * block_width;
          tab_len = 1 << block_width;
        }

      /* Allocate and initialize a table. */
      table->tabsx[i] = (mpz_t *)malloc(tab_len * sizeof(mpz_t));
      tx = table->tabsx[i];
      table->tabsy[i] = (mpz_t *)malloc(tab_len * sizeof(mpz_t));
      ty = table->tabsy[i];

      /* Initialize mpz_t's. */
      for (j = 0; j < tab_len; j++)
        {
          mpz_init(tx[j]);
          mpz_init(ty[j]);
        }
    }
}
