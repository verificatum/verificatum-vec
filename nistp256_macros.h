
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

#include "undefine_macros.h"

#define POSTFIX _nistp256_inner
#define TAB_POSTFIX _nistp256_inner

#define FIELD_ELEMENT felem
#define FIELD_ELEMENT_INIT(x)
#define FIELD_ELEMENT_CLEAR(x)
#define FIELD_ELEMENT_UNIT(x, y, z) \
  memset(x, 0, sizeof(felem));      \
  memset(y, 0, sizeof(felem));      \
  memset(z, 0, sizeof(felem))
#define FIELD_ELEMENT_SET(rx, ry, rz, x, y, z) \
  felem_assign(rx, x);                         \
  felem_assign(ry, y);                         \
  felem_assign(rz, z)

#define FIELD_ELEMENT_VAR smallfelem
#define FIELD_ELEMENT_VAR_INIT(x)
#define FIELD_ELEMENT_VAR_CLEAR(x)
#define FIELD_ELEMENT_VAR_UNIT(x, y, z) \
  memset(x, 0, sizeof(smallfelem));     \
  memset(y, 0, sizeof(smallfelem));     \
  memset(z, 0, sizeof(smallfelem))
#define FIELD_ELEMENT_VAR_SET(rx, ry, rz, x, y, z) \
  smallfelem_assign(rx, x);                        \
  smallfelem_assign(ry, y);                        \
  smallfelem_assign(rz, z)

#define FIELD_ELEMENT_CONTRACT(rx, ry, rz, x, y, z) \
  felem_contract(rx, x);                            \
  felem_contract(ry, y);                            \
  felem_contract(rz, z);


#define SCRATCH(scratch)
#define SCRATCH_INIT(scratch)
#define SCRATCH_CLEAR(scratch)

#define ARRAY_MALLOC_INIT(len) \
  (FIELD_ELEMENT_VAR *)malloc(len * sizeof(FIELD_ELEMENT_VAR))

#define ARRAY_CLEAR_FREE(array, len) free(array)

#define JDBL(scratch, rx, ry, rz, curve, x, y, z) \
  point_double(rx, ry, rz, x, y, z)

#define JDBL_VAR(scratch, rx, ry, rz, curve, x, y, z) \
  point_double_small(rx, ry, rz, x, y, z)             \

#define JADD(scratch, rx, ry, rz, curve, x1, y1, z1, x2, y2, z2) \
  point_add(rx, ry, rz, x1, y1, z1, x2, y2, z2);

#define JADD_VAR(scratch, rx, ry, rz, curve, x1, y1, z1, x2, y2, z2)  \
  point_add_small(rx, ry, rz, x1, y1, z1, x2, y2, z2)

#define CURVE vec_curve
