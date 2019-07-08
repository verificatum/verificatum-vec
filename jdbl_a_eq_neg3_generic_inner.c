
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
#define alpha scratch->t3
#define beta scratch->t4
#define gamma scratch->t5
#define delta scratch->t6

#define modulus curve->modulus
#define a curve->a


/* 2001 Bernstein Jacobi coordinates. Special case a = -3. */

void
vec_jdbl_a_eq_neg3_generic_inner(vec_scratch_mpz_t scratch,
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

  /* delta = Z1^2 */
  mpz_mul(delta, Z1, Z1);
  mpz_mod(delta, delta, modulus);

  /* gamma =Y1^2 */
  mpz_mul(gamma, Y1, Y1);
  mpz_mod(gamma, gamma, modulus);

  /* beta = X1*gamma */
  mpz_mul(beta, X1, gamma);
  mpz_mod(beta, beta, modulus);

  /* alpha = 3*(X1-delta)*(X1+delta) */
  mpz_sub(t1, X1, delta);
  mpz_add(t2, X1, delta);
  mpz_mul_si(t1, t1, 3);
  mpz_mul(alpha, t1, t2);
  mpz_mod(alpha, alpha, modulus);

  /* X3 = alpha^2-8*beta */
  mpz_mul(t1, alpha, alpha);
  mpz_mul_si(t2, beta, 8);
  mpz_sub(X3, t1, t2);
  mpz_mod(X3, X3, modulus);

  /* Z3 = (Y1+Z1)^2-gamma-delta */
  mpz_add(t1, Y1, Z1);
  mpz_mul(t1, t1, t1);
  mpz_sub(t1, t1, gamma);
  mpz_sub(t1, t1, delta);
  mpz_mod(Z3, t1, modulus);

  /* Y3 = alpha*(4*beta-X3)-8*gamma^2 */
  mpz_mul_si(t1, beta, 4);
  mpz_sub(t1, t1, X3);
  mpz_mul(t1, t1, alpha);

  mpz_mul(t2, gamma, gamma);
  mpz_mul_si(t2, t2, 8);

  mpz_sub(Y3, t1, t2);
  mpz_mod(Y3, Y3, modulus);
}
