
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

#ifndef JSMUL_H_TEMPLATE_H
#define JSMUL_H_TEMPLATE_H

#include <gmp.h>
#include "vec.h"
#include "templates.h"

/**
 * Simultaneous multiplication table.
 */
typedef struct
{
  vec_curve *curve;           /**< Underlying curve. */
  size_t len;                 /**< Total number of bases/scalars. */
  size_t block_width;         /**< Number of bases/scalars in each block. */
  size_t tabs_len;            /**< Number of blocks. */
  FIELD_ELEMENT_VAR **tabsx;  /**< Table of tables, one sub-table for each
                                block. */
  FIELD_ELEMENT_VAR **tabsy;  /**< Table of tables, one sub-table for each
                                 block. */
  FIELD_ELEMENT_VAR **tabsz;  /**< Table of tables, one sub-table for each
                                 block. */

} FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX)[1]; /* Magic references. */

void
FUNCTION_NAME(vec_jsmul_init, POSTFIX)
     (FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX) table,
      CURVE *curve,
      size_t len,
      size_t block_width);

int
FUNCTION_NAME(vec_jsmul_clear, POSTFIX)
     (FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX) table);

void
FUNCTION_NAME(vec_jsmul_precomp, POSTFIX)
     (FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX) table,
      CURVE *curve,
      FIELD_ELEMENT_VAR *basesx, FIELD_ELEMENT_VAR *basesy,
      FIELD_ELEMENT_VAR *basesz);

void
FUNCTION_NAME(vec_jsmul_table, POSTFIX)
     (FIELD_ELEMENT_VAR ropx, FIELD_ELEMENT_VAR ropy, FIELD_ELEMENT_VAR ropz,
      CURVE *curve,
      FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX) table,
      mpz_t *scalars,
      size_t max_scalar_bitlen);

void
FUNCTION_NAME(vec_jsmul_block_batch, POSTFIX)
     (FIELD_ELEMENT ropx, FIELD_ELEMENT ropy, FIELD_ELEMENT ropz,
      CURVE *curve,
      FIELD_ELEMENT_VAR *basesx, FIELD_ELEMENT_VAR *basesy,
      FIELD_ELEMENT_VAR *basesz,
      mpz_t *scalars,
      size_t len,
      size_t block_width, size_t batch_len,
      size_t max_scalar_bitlen);

#endif /* JSMUL_H_TEMPLATE_H */
