
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

#include <stdlib.h>
#include "templates.h"

#include "jfmul_h_template.h"
#include "jsmul_h_template.h"

void
FUNCTION_NAME(vec_jfmul_init, POSTFIX)
     (FUNCTION_NAME(vec_jfmul_tab, TAB_POSTFIX) *table,
      CURVE *curve,
      size_t len)
{
  size_t bit_length = mpz_sizeinbase(curve->n, 2);

  int block_width = vec_fmul_block_width(bit_length, len);

  FUNCTION_NAME(vec_jsmul_init, POSTFIX)(table->tab,
                                           curve,
                                           block_width,
                                           block_width);
  table->slice_bit_len =
    (((int)bit_length) + (block_width - 1)) / block_width;
}

void
FUNCTION_NAME(vec_jfmul_clear_free, POSTFIX)
     (FUNCTION_NAME(vec_jfmul_tab, TAB_POSTFIX) *table)
{
  FUNCTION_NAME(vec_jsmul_clear, POSTFIX)(table->tab);
  free(table);
}

void
FUNCTION_NAME(vec_jfmul_prcmp, POSTFIX)
     (CURVE *curve,
      FUNCTION_NAME(vec_jfmul_tab, TAB_POSTFIX) *table,
      FIELD_ELEMENT_VAR x, FIELD_ELEMENT_VAR y, FIELD_ELEMENT_VAR z)
{
  size_t i;
  size_t j;
  size_t bw = table->tab->block_width;

  FIELD_ELEMENT_VAR *basesx;
  FIELD_ELEMENT_VAR *basesy;
  FIELD_ELEMENT_VAR *basesz;

  SCRATCH(scratch);

  SCRATCH_INIT(scratch);

  basesx = ARRAY_MALLOC_INIT(bw);
  basesy = ARRAY_MALLOC_INIT(bw);
  basesz = ARRAY_MALLOC_INIT(bw);

  FIELD_ELEMENT_VAR_SET(basesx[0], basesy[0], basesz[0], x, y, z);

  for (i = 1; i < bw; i++)
    {

      FIELD_ELEMENT_VAR_SET(basesx[i], basesy[i], basesz[i],
                            basesx[i - 1], basesy[i - 1], basesz[i - 1]);

      for (j = 0; j < table->slice_bit_len; j++)
        {
          JDBL_VAR(scratch,
                   basesx[i], basesy[i], basesz[i],
                   curve,
                   basesx[i], basesy[i], basesz[i]);
        }
    }

  FUNCTION_NAME(vec_jsmul_precomp, POSTFIX)(table->tab,
                                              curve,
                                              basesx, basesy, basesz);

  ARRAY_CLEAR_FREE(basesx, bw);
  ARRAY_CLEAR_FREE(basesy, bw);
  ARRAY_CLEAR_FREE(basesz, bw);

  SCRATCH_CLEAR(scratch);
}

void
FUNCTION_NAME(vec_jfmul_cmp, POSTFIX)
     (FIELD_ELEMENT_VAR ropx, FIELD_ELEMENT_VAR ropy, FIELD_ELEMENT_VAR ropz,
      CURVE *curve, FUNCTION_NAME(vec_jfmul_tab, TAB_POSTFIX) *table,
      mpz_t scalar)
{
  size_t i;
  size_t bw = table->tab->block_width;
  mpz_t *scalars;
  mpz_t tmp;

  mpz_init(tmp);

  scalars = vec_array_alloc_init(bw);

  mpz_set(tmp, scalar);

  for (i = 0; i < bw; i++)
    {
      mpz_tdiv_r_2exp(scalars[i], tmp, table->slice_bit_len);
      mpz_tdiv_q_2exp(tmp, tmp, table->slice_bit_len);
    }

  FUNCTION_NAME(vec_jsmul_table, POSTFIX)(ropx, ropy, ropz,
                                            curve,
                                            table->tab,
                                            scalars,
                                            table->slice_bit_len);
  vec_array_clear_free(scalars, bw);

  mpz_clear(tmp);
}
