
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

/* Naive version of multiplication. Only used during development. */

/* LCOV_EXCL_START */
#ifdef FIELD_ELEMENT

#include <stdio.h>
#include "templates.h"

void
FUNCTION_NAME(vec_jmul, POSTFIX)
     (FIELD_ELEMENT RX, FIELD_ELEMENT RY, FIELD_ELEMENT RZ,
      CURVE *curve,
      FIELD_ELEMENT_VAR X, FIELD_ELEMENT_VAR Y, FIELD_ELEMENT_VAR Z,
      mpz_t scalar)
{
  int i;
  int bitLength;

  SCRATCH(scratch);

  VEC_UNUSED(curve);

  SCRATCH_INIT(scratch);

  FIELD_ELEMENT_UNIT(RX, RY, RZ);

  /* Determine bit length. */
  bitLength = (int)mpz_sizeinbase(scalar, 2);

  for (i = bitLength; i >= 0; i--)
    {

      /* Double. */
      JDBL(scratch,
           RX, RY, RZ,
           curve,
           RX, RY, RZ);

      /* Add. */
      if (mpz_tstbit(scalar, i))
        {

          JADD(scratch,
               RX, RY, RZ,
               curve,
               RX, RY, RZ,
               X, Y, Z);
        }
    }

  SCRATCH_CLEAR(scratch);
}

#endif
/* LCOV_EXCL_STOP */
