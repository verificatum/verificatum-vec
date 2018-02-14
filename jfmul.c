
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

#include "generic_macros.h"
#include "jfmul_template.h"

size_t
vec_jfmul_ptrsize()
{
  return sizeof(vec_jfmul_tab*);
}

vec_jfmul_tab_ptr
vec_jfmul_precomp(vec_curve *curve,
                    mpz_t X, mpz_t Y, mpz_t Z,
                    size_t len)
{
  vec_jfmul_tab_ptr ptr;

  ptr.generic = (vec_jfmul_tab*)malloc(sizeof(vec_jfmul_tab));

  vec_jfmul_init(ptr.generic, curve, len);

  vec_jfmul_prcmp(curve, ptr.generic, X, Y, Z);

  return ptr;
}

void
vec_jfmul(mpz_t RX, mpz_t RY, mpz_t RZ,
            vec_curve *curve, vec_jfmul_tab_ptr table_ptr,
            mpz_t scalar)
{

  vec_jfmul_cmp(RX, RY, RZ,
                  curve, table_ptr.generic,
                  scalar);
}

void
vec_jfmul_free(vec_jfmul_tab_ptr table_ptr)
{
  vec_jfmul_clear_free(table_ptr.generic);
}
