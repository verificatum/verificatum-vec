
/*
 * Copyright 2008-2018 Douglas Wikstrom
 *
 * This file is part of Verificatum Elliptic Curve library (VEC).
 *
 * VEC is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * VEC is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with VEC. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <gmp.h>
#include "vec.h"

#define t1 scratch->t1
#define t2 scratch->t2
#define t3 scratch->t3
#define U1 scratch->t4
#define U2 scratch->t5
#define S1 scratch->t6
#define S2 scratch->t7
#define H scratch->t8
#define r scratch->t9

#define modulus curve->modulus
#define a curve->a

/* 1998 Cohen/Miyaji/Ono Jacobi coordinates with cached powers of
   Z2. */

void
vec_jadd_generic_inner(vec_scratch_mpz_t scratch,
                       mpz_t X3, mpz_t Y3, mpz_t Z3,
                       vec_curve *curve,
                       mpz_t X1, mpz_t Y1, mpz_t Z1,
                       mpz_t X2, mpz_t Y2, mpz_t Z2)
{

  /* P1 is point at infinity. */
  if (mpz_cmp_si(Z1, 0) == 0)
    {

      /* P2 is also point at infinity. */
      if (mpz_cmp_si(Z2, 0) == 0)
        {
          mpz_set_si(X3, 0);
          mpz_set_si(Y3, 1);
          mpz_set_si(Z3, 0);
          return;
        }
      /* P1 is point at infinity and P2 is not. */
      else
        {
          mpz_set(X3, X2);
          mpz_set(Y3, Y2);
          mpz_set(Z3, Z2);
          return;
        }
    }

  /* P2 is point at infinity and P1 is not. */
  else if (mpz_cmp_si(Z2, 0) == 0)
    {

      mpz_set(X3, X1);
      mpz_set(Y3, Y1);
      mpz_set(Z3, Z1);
      return;
    }

  /* Compute powers of Z2. */
  mpz_mul(t1, Z2, Z2);           /* t1 = Z2^2 */
  mpz_mod(t1, t1, modulus);
  mpz_mul(S2, t1, Z2);           /* S2 = Z2^3 */
  mpz_mod(S2, S2, modulus);

  /* Compute powers of Z1 */
  mpz_mul(t2, Z1, Z1);           /* t2 = Z1^2 */
  mpz_mod(t2, t2, modulus);
  mpz_mul(t3, t2, Z1);           /* t3 = Z1^3 */
  mpz_mod(t3, t3, modulus);

  /* U1:=X1*Z2^2 */
  mpz_mul(U1, X1, t1);
  mpz_mod(U1, U1, modulus);

  /* U2:=X2*Z1^2 */
  mpz_mul(U2, X2, t2);

  /* S1:=Y1*Z2^3 */
  mpz_mul(S1, Y1, S2);
  mpz_mod(S1, S1, modulus);

  /* S2:=Y2*Z1^3 */
  mpz_mul(S2, Y2, t3);

  /* H:=U2-U1 */
  mpz_sub(H, U2, U1);
  mpz_mod(H, H, modulus);

  /* r:=S2-S1 */
  mpz_sub(r, S2, S1);
  mpz_mod(r, r, modulus);

  if (mpz_cmp_si(H, 0) == 0)
    {

      if (mpz_cmp_si(r, 0) != 0)
        {
          mpz_set_si(X3, 0);
          mpz_set_si(Y3, 1);
          mpz_set_si(Z3, 0);
          return;
        }
      else
        {

          curve->jdbl(scratch,
                      X3, Y3, Z3,
                      curve,
                      X1, Y1, Z1);
          return;
        }
    }

  /* Compute square of r */
  mpz_mul(t1, r, r);          /* t1 = r^2 */
  mpz_mod(t1, t1, modulus);

  /* Compute powers of H */
  mpz_mul(t2, H, H);          /* t2 = H^2 */
  mpz_mod(t2, t2, modulus);
  mpz_mul(t3, t2, H);         /* t3 = H^3 */
  mpz_mod(t3, t3, modulus);


  /* X3:=-H^3-2*U1*H^2+r^2 */
  mpz_sub(X3, t1, t3);        /* X3 = r^2 - H^3 */

  mpz_mul(t1, U1, t2);        /* t1 = 2*U1*H^2 */
  mpz_mul_si(t1, t1, 2);
  mpz_mod(t1, t1, modulus);

  mpz_sub(X3, X3, t1);
  mpz_mod(X3, X3, modulus);

  /* Y3:=-S1*H^3+r*(U1*H^2-X3) */
  mpz_mul(t1, U1, t2);        /* t1 = r*(U1*H^2-X3) */
  mpz_mod(t1, t1, modulus);
  mpz_sub(t1, t1, X3);
  mpz_mul(t1, r, t1);
  mpz_mod(t1, t1, modulus);

  mpz_mul(t2, S1, t3);        /* t2 = S1*H^3 */
  mpz_mod(t2, t2, modulus);

  mpz_sub(Y3, t1, t2);
  mpz_mod(Y3, Y3, modulus);

  /* Z3:=Z1*Z2*H */
  mpz_mul(Z3, Z1, Z2);
  mpz_mod(Z3, Z3, modulus);
  mpz_mul(Z3, Z3, H);
  mpz_mod(Z3, Z3, modulus);
}
