
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
vec_jsmul_aff(mpz_t ropx, mpz_t ropy,
              vec_curve *curve,
              mpz_t *basesx, mpz_t *basesy,
              mpz_t *exponents,
              size_t len)
{
  size_t i;

  mpz_t ropz;
  mpz_t *basesz = vec_array_alloc_init(len);

  mpz_init(ropz);

  for (i = 0; i < len; i++) {
    vec_affj(basesx[i], basesy[i], basesz[i]);
  }

  curve->jsmul(ropx, ropy, ropz,
               curve,
               basesx, basesy, basesz,
               exponents,
               len);

  vec_jaff(ropx, ropy, ropz, curve);

  mpz_clear(ropz);

  vec_array_clear_free(basesz, len);
}
