
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

void
vec_jfmul_aff(mpz_t rx, mpz_t ry,
                vec_curve *curve,
                vec_jfmul_tab_ptr table_ptr,
                mpz_t scalar)
{

  mpz_t RZ;

  mpz_init(RZ);

  curve->jfmul(rx, ry, RZ,
               curve, table_ptr,
               scalar);

  vec_jaff(rx, ry, RZ, curve);

  mpz_clear(RZ);
}
