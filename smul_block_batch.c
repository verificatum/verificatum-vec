
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
