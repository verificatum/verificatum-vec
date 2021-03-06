
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

#define t1 scratch->t1
#define t2 scratch->t2
#define s scratch->t3

#define modulus curve->modulus
#define a curve->a
#define b curve->b

void
vec_add(vec_scratch_mpz_t scratch,
        mpz_t rx, mpz_t ry,
        vec_curve *curve,
        mpz_t x1, mpz_t y1,
        mpz_t x2, mpz_t y2)
{

  /* If the first point is the unit element, then we copy the other
     point. */
  if (mpz_cmp_si(x1, -1) == 0)
    {
      mpz_set(rx, x2);
      mpz_set(ry, y2);
      return;
    }

  /* If the second point is the unit element, then we copy the other
     point. */
  if (mpz_cmp_si(x2, -1) == 0)
    {
      mpz_set(rx, x1);
      mpz_set(ry, y1);
      return;
    }

  /* If the second point is inverse of the first point, then we return
     the unit point. */
  if (mpz_cmp(x1, x2) == 0)
    {

      mpz_add(t1, y1, y2);

      if (mpz_cmp(t1, modulus) == 0)
        {
          mpz_set_si(rx, -1);
          mpz_set_si(ry, -1);
          return;
        }
    }

  /* If the first and second points are identical, then we double
     it. */
  if (mpz_cmp(x1, x2) == 0 && mpz_cmp(y1, y2) == 0)
    {
      vec_dbl(scratch, rx, ry, curve, x1, y1);
      return;
    }

  /* s = (y1 - y2) / (x1 - x2) */
  mpz_sub(t1, y1, y2);
  mpz_sub(t2, x1, x2);
  mpz_invert(t2, t2, modulus);
  mpz_mul(s, t1, t2);
  mpz_mod(s, s, modulus);

  /* rx = s^2 - (x1 + x2) */
  mpz_mul(t1, s, s);
  mpz_sub(t1, t1, x1);
  mpz_sub(t1, t1, x2);

  /* ry = s(x1 - rx) - y1 */
  mpz_sub(t2, x1, t1);
  mpz_mul(t2, s, t2);
  mpz_sub(t2, t2, y1);

  /* We assign the destination parameters in the end to allow them to
     be identical to the inputs. */
  mpz_mod(rx, t1, modulus);
  mpz_mod(ry, t2, modulus);

}
