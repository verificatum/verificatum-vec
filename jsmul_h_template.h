
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
