
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
