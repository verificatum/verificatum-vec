
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
