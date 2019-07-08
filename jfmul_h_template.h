
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

#ifndef JFMUL_H_TEMPLATE_H
#define JFMUL_H_TEMPLATE_H

#include <gmp.h>
#include "vec.h"
#include "templates.h"
#include "jsmul_h_template.h"

struct FUNCTION_NAME(_vec_jfmul_tab, TAB_POSTFIX)
{

  FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX) tab; /**< Underlying simultaneous
                                                multiplication table. */
  size_t slice_bit_len;                      /**< Bit length of each slice. */

};
typedef struct FUNCTION_NAME(_vec_jfmul_tab, TAB_POSTFIX)
FUNCTION_NAME(vec_jfmul_tab, TAB_POSTFIX); /* Magic references. */

void
FUNCTION_NAME(vec_jfmul_init, POSTFIX)
     (FUNCTION_NAME(vec_jfmul_tab, TAB_POSTFIX) *table,
      CURVE *curve,
      size_t len);

void
FUNCTION_NAME(vec_jfmul_clear_free, POSTFIX)
     (FUNCTION_NAME(vec_jfmul_tab, TAB_POSTFIX) *table);

void
FUNCTION_NAME(vec_jfmul_prcmp, POSTFIX)
     (CURVE *curve,
      FUNCTION_NAME(vec_jfmul_tab, TAB_POSTFIX) *table,
      FIELD_ELEMENT_VAR x, FIELD_ELEMENT_VAR y, FIELD_ELEMENT_VAR z);

void
FUNCTION_NAME(vec_jfmul_cmp, POSTFIX)
     (FIELD_ELEMENT_VAR ropx, FIELD_ELEMENT_VAR ropy, FIELD_ELEMENT_VAR ropz,
      CURVE *curve, FUNCTION_NAME(vec_jfmul_tab, TAB_POSTFIX) *table,
      mpz_t scalar);


#endif /* JFMUL_H_TEMPLATE_H */
