
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
