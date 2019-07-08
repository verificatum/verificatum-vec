
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
#define t3 scratch->t3
#define S scratch->t4
#define M scratch->t5
#define T scratch->t6
#define modulus curve->modulus
#define a curve->a

/* 1998 Cohen/Miyaji/Ono Jacobi coordinates. */

void
vec_jdbl_generic_inner(vec_scratch_mpz_t scratch,
                       mpz_t X3, mpz_t Y3, mpz_t Z3,
                       vec_curve *curve,
                       mpz_t X1, mpz_t Y1, mpz_t Z1)
{

  /* (X1, Y1, Z1) is point at infinity or point which is its own
     inverse. */
  if (mpz_cmp_ui(Z1, 0) == 0 || mpz_cmp_ui(Y1, 0) == 0)
    {
      mpz_set_ui(X3, 0);
      mpz_set_ui(Y3, 1);
      mpz_set_ui(Z3, 0);
      return;
    }

  /* S = 4*X1*Y1^2 */
  mpz_mul(S, Y1, Y1);
  mpz_mod(S, S, modulus);
  mpz_mul(S, S, X1);
  mpz_mul_si(S, S, 4);
  mpz_mod(S, S, modulus);


  /* Z1 squared */
  mpz_mul(t2, Z1, Z1);          /* t2 = Z1^2 */
  mpz_mod(t2, t2, modulus);

  /* M = 3*X1^2+a*Z1^4 */
  mpz_mul(t1, X1, X1);          /* t1 = 3*X1^2 */
  mpz_mod(t1, t1, modulus);
  mpz_mul_si(t1, t1, 3);
  mpz_mod(t1, t1, modulus);

  mpz_mul(t3, t2, t2);          /* t3 = a*Z1^4 */
  mpz_mod(t3, t3, modulus);
  mpz_mul(t3, t3, a);
  mpz_mod(t3, t3, modulus);

  mpz_add(M, t1, t3);
  mpz_mod(M, M, modulus);

  /* T = M^2-2*S */
  mpz_mul(T, M, M);
  mpz_mul_si(t2, S, 2);
  mpz_sub(T, T, t2);
  mpz_mod(T, T, modulus);

  /* X3 = T */
  mpz_set(X3, T);

  /* Y3 = -8*Y1^4+M*(S-T) */
  mpz_sub(t1, S, T);            /* t1 = M*(S-T) */
  mpz_mul(t1, t1, M);
  mpz_mod(t1, t1, modulus);

  mpz_mul(t2, Y1, Y1);          /* t2 = 8*Y1^4 */
  mpz_mod(t2, t2, modulus);
  mpz_mul(t2, t2, t2);
  mpz_mod(t2, t2, modulus);
  mpz_mul_si(t2, t2, 8);
  mpz_mod(t2, t2, modulus);

  mpz_sub(t1, t1, t2);

  /* Z3 = 2*Y1*Z1 */
  mpz_mul(t2, Y1, Z1);
  mpz_mul_si(t2, t2, 2);

  mpz_mod(Y3, t1, modulus);
  mpz_mod(Z3, t2, modulus);

}
