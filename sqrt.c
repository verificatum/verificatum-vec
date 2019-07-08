
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

void
vec_sqrt(mpz_t res, mpz_t a, mpz_t p) {

  int s;
  int t;

  mpz_t v;
  mpz_t k;
  mpz_t r;
  mpz_t n;
  mpz_t z;
  mpz_t c;

  mpz_init(v);
  mpz_init(k);
  mpz_init(r);
  mpz_init(n);
  mpz_init(z);
  mpz_init(c);

  /* LCOV_EXCL_START */
  /* Square root of zero is zero. */
  if (mpz_cmp_si(a, 0) == 0) {
    mpz_set_si(res, 0);
    goto CLEAR;
  }
  /* LCOV_EXCL_STOP */

  /* If p = 3 mod 4, then computing a square root is trivial. */
  if (mpz_tstbit(p, 0) == 1 && mpz_tstbit(p, 1) == 1) {

    mpz_add_ui(v, p, 1);           /* v = (p + 1)/4 */
    mpz_tdiv_q_2exp(v, v, 2);
    mpz_powm(res, a, v, p);        /* res = a^((p+1)/4) */

    goto CLEAR;
  }


  /* Compute k and s, where p = 2^s(2k+1) + 1 */
  s = 0;
  mpz_sub_ui(k, p, 1);
  while (mpz_tstbit(k, 0) == 0) {
    s++;
    mpz_tdiv_q_2exp(k, k, 1);
  }
  mpz_sub_ui(k, k, 1);
  mpz_tdiv_q_2exp(k, k, 1);

  mpz_powm(r, a, k, p);  /* r = a^k mod p */

  mpz_mul(n, r, r);      /* n = r^2 * a mod p */
  mpz_mod(n, n, p);
  mpz_mul(n, n, a);
  mpz_mod(n, n, p);

  mpz_mul(r, r, a);      /* r = r * a mod p*/
  mpz_mod(r, r, p);

  if (mpz_cmp_si(n, 1) == 0) {
    mpz_set(res, r);
    goto CLEAR;
  }

  mpz_set_si(z, 2);      /* z = quadratic non-residue */
  while (mpz_legendre(z, p) == 1) {
    mpz_add_ui(z, z, 1);
  }

  mpz_set(v, k);         /* v = 2k + 1 */
  mpz_mul_si(v, v, 2);
  mpz_add_ui(v, v, 1);

  /* c = z^v mod p */
  mpz_powm(c, z, v, p);

  /* Iterate */
  while (mpz_cmp_si(n, 1) > 0) {

    mpz_set(k, n);
    t = s;
    s = 0;

    while (mpz_cmp_si(k, 1) != 0) {

      mpz_mul(k, k, k);  /* k = k^2 mod p */
      mpz_mod(k, k, p);

      s++;
    }

    t -= s;

    mpz_set_si(v, 1);      /* v = 2^(t-1) */
    mpz_mul_2exp(v, v, t - 1);

    mpz_powm(c, c, v, p);  /* c = c^v mod p */
    mpz_mul(r, r, c);      /* r = rc mod p  */
    mpz_mod(r, r, p);
    mpz_mul(c, c, c);      /* c = c^2 mod p */
    mpz_mod(c, c, p);
    mpz_mul(n, n, c);      /* n = nc mod p  */
    mpz_mod(n, n, p);
  }

  mpz_mod(res, r, p);


CLEAR:

  mpz_clear(c);
  mpz_clear(z);
  mpz_clear(n);
  mpz_clear(r);
  mpz_clear(k);
  mpz_clear(v);

}
