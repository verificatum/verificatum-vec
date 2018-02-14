
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
