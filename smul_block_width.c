
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

/*
 * Computes a "theoretical" optimal block width for a given scalar
 * length.
 */
int
vec_smul_block_width(int scalars_bitlen, int batch_len) {

  /* This computes the theoretical optimum. */
  int width = 1;
  double cost = 1.5 * scalars_bitlen;
  double dbl;
  double add_precomp;
  double add;
  double old_cost;
  int width_exp;

  do {

    old_cost = cost;

    width++;
    width_exp = 1 << width;

    /* Amortized cost for doublings. */
    dbl = ((double)scalars_bitlen) / batch_len;

    /* Amortized cost for precomputing for a block. */
    add_precomp = ((double)width_exp) / width;

    /* Amortized cost for adding, given precomputation. */
    add = (((double)(scalars_bitlen)) / width) * (1 - 1.0 / width_exp);

    /* Total amortized cost. */
    cost = dbl + add_precomp + add;

  } while (cost < old_cost);

  width--;

  return width;
}
