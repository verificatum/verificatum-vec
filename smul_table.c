
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
#include "gmp.h"
#include "vec.h"

/*
 * Returns the index'th bit of each of the first block_width integers
 * in the array. The least significant bit in the output is the bit
 * extracted from the first integer in the input array.
 */
static int
getbits(mpz_t *op, int index, size_t block_width)
{
  int i;
  int bits = 0;

  for (i = block_width - 1; i >= 0; i--)
    {
      bits <<= 1;
      if (mpz_tstbit(op[i], index))
        {
          bits |= 1;
        }
    }
  return bits;
}

void
vec_smul_table(mpz_t ropx, mpz_t ropy,
               vec_smul_tab table,
               mpz_t *scalars,
               size_t max_scalar_bitlen)
{
  size_t i;
  int index;
  int mask;

  vec_scratch_mpz_t scratch;
  mpz_t *exps;

  size_t len = table->len;
  size_t tabs_len = table->tabs_len;
  size_t block_width = table->block_width;
  size_t last_block_width = len - (tabs_len - 1) * block_width;
  mpz_t **tabsx = table->tabsx;
  mpz_t **tabsy = table->tabsy;

  vec_scratch_init_mpz_t(scratch);

  /* Initialize result variable. */
  mpz_set_si(ropx, -1);
  mpz_set_si(ropy, -1);

  /* Execute simultaneous double-and-add. */
  for (index = max_scalar_bitlen - 1; index >= 0; index--)
    {

      /* Double ... */
      vec_dbl(scratch,
              ropx, ropy,
              table->curve,
              ropx, ropy);

      /* ... and add. */
      i = 0;
      exps = scalars;
      while (i < tabs_len)
        {
          if (i == tabs_len - 1)
            {
              mask = getbits(exps, index, last_block_width);
            }
          else
            {
              mask = getbits(exps, index, block_width);
            }

          vec_add(scratch,
                  ropx, ropy,
                  table->curve,
                  ropx, ropy,
                  tabsx[i][mask], tabsy[i][mask]);
          i++;
          exps += block_width;
        }
    }

  vec_scratch_clear_mpz_t(scratch);
}
