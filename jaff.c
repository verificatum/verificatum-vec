
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

#include <gmp.h>
#include "vec.h"

void
vec_jaff(mpz_t X, mpz_t Y, mpz_t Z, vec_curve *curve) {

  mpz_t ZZ;

  if (mpz_cmp_si(Z, 0) == 0)
    {
      mpz_set_si(X, -1);
      mpz_set_si(Y, -1);
    }
  else
    {

      mpz_init(ZZ);

      mpz_invert(Z, Z, curve->modulus);

      mpz_mul(ZZ, Z, Z);
      mpz_mod(ZZ, ZZ, curve->modulus);

      mpz_mul(X, X, ZZ);
      mpz_mod(X, X, curve->modulus);

      mpz_mul(Z, ZZ, Z);
      mpz_mod(Z, Z, curve->modulus);

      mpz_mul(Y, Y, Z);
      mpz_mod(Y, Y, curve->modulus);

      mpz_clear(ZZ);
    }
}
