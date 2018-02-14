
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
vec_smul_clear(vec_smul_tab table)
{
  size_t i, j;                             /* Index variables. */
  size_t tabs_len = table->tabs_len;       /* Number of sub tables. */
  size_t block_width = table->block_width; /* Width of each sub table. */
  size_t tab_len = 1 << block_width;       /* Size of each sub table. */
  mpz_t *tx;                               /* Temporary table variable. */
  mpz_t *ty;                               /* Temporary table variable. */


  for (i = 0; i < tabs_len; i++)
    {

      /* Last block may have smaller width, but the width is never
         zero. */
      if (i == tabs_len - 1)
        {
          block_width = table->len - (tabs_len - 1) * block_width;
          tab_len = 1 << block_width;
        }

      /* Deallocate all integers in table. */
      tx = table->tabsx[i];
      ty = table->tabsy[i];
      for (j = 0; j < tab_len; j++)
        {
          mpz_clear(tx[j]);
          mpz_clear(ty[j]);
        }

      /* Deallocate table. */
      free(tx);
      free(ty);
    }

  /* Deallocate table of tables. */
  free(table->tabsx);
  free(table->tabsy);
}
