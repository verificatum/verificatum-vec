
/* Copyright 2008-2019 Douglas Wikstrom
 *
 * This file is part of Verificatum Elliptic Curve library (VEC).
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
