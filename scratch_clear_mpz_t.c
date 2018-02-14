
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

#include <gmp.h>

#include "vec.h"

void
vec_scratch_clear_mpz_t(vec_scratch_mpz_t scratch)
{
  mpz_clear(scratch->t1);
  mpz_clear(scratch->t2);
  mpz_clear(scratch->t3);
  mpz_clear(scratch->t4);
  mpz_clear(scratch->t5);
  mpz_clear(scratch->t6);
  mpz_clear(scratch->t7);
  mpz_clear(scratch->t8);
  mpz_clear(scratch->t9);
}
