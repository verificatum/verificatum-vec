
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

vec_curve *
vec_curve_get(mpz_t modulus, mpz_t a, mpz_t b,
                mpz_t gx, mpz_t gy, mpz_t n)
{
  vec_curve *curve = vec_curve_alloc();

  curve->name = NULL;

  mpz_set(curve->modulus, modulus);
  mpz_set(curve->a, a);
  mpz_set(curve->b, b);

  mpz_set(curve->gx, gx);
  mpz_set(curve->gy, gy);
  mpz_set(curve->n, n);

  curve->jdbl = vec_jdbl;
  curve->jadd = vec_jadd;
  curve->jmul = vec_jmul;
}
