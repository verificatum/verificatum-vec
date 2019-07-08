
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
vec_jdbl_aff(vec_scratch_mpz_t scratch,
             mpz_t rx, mpz_t ry,
             vec_curve *curve,
             mpz_t x, mpz_t y) {
  mpz_t X1;
  mpz_t Y1;
  mpz_t Z1;
  mpz_t Z3;

  mpz_init(X1);
  mpz_init(Y1);
  mpz_init(Z1);
  mpz_init(Z3);

  mpz_set(X1, x);
  mpz_set(Y1, y);

  vec_affj(X1, Y1, Z1);

  curve->jdbl(scratch,
              rx, ry, Z3,
              curve,
              X1, Y1, Z1);

  vec_jaff(rx, ry, Z3, curve);

  mpz_clear(Z3);
  mpz_clear(Z1);
  mpz_clear(Y1);
  mpz_clear(X1);
}
