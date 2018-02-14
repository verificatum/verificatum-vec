
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

#include "undefine_macros.h"

#define POSTFIX _nistp521_inner
#define TAB_POSTFIX _nistp521_inner

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

#define FIELD_ELEMENT_VAR felem
#define FIELD_ELEMENT_VAR_INIT(x)
#define FIELD_ELEMENT_VAR_CLEAR(x)
#define FIELD_ELEMENT_VAR_UNIT(x, y, z) \
  memset(x, 0, sizeof(felem));          \
  memset(y, 0, sizeof(felem));          \
  memset(z, 0, sizeof(felem))
#define FIELD_ELEMENT_VAR_SET(rx, ry, rz, x, y, z) \
  felem_assign(rx, x);                             \
  felem_assign(ry, y);                             \
  felem_assign(rz, z)

#define FIELD_ELEMENT_CONTRACT(rx, ry, rz, x, y, z) \
  felem_assign(rx, x);                              \
  felem_assign(ry, y);                              \
  felem_assign(rz, z)

#define SCRATCH(scratch)
#define SCRATCH_INIT(scratch)
#define SCRATCH_CLEAR(scratch)

#define ARRAY_MALLOC_INIT(len) \
  (FIELD_ELEMENT *)malloc(len * sizeof(FIELD_ELEMENT))

#define ARRAY_CLEAR_FREE(array, len) free(array)

#define JDBL(scratch, rx, ry, rz, curve, x, y, z) \
  point_double(rx, ry, rz, x, y, z)

#define JDBL_VAR(scratch, rx, ry, rz, curve, x, y, z) \
  point_double(rx, ry, rz, x, y, z)

#define JADD(scratch, rx, ry, rz, curve, x1, y1, z1, x2, y2, z2) \
  point_add(rx, ry, rz, x1, y1, z1, x2, y2, z2)

#define JADD_VAR(scratch, rx, ry, rz, curve, x1, y1, z1, x2, y2, z2) \
  point_add(rx, ry, rz, x1, y1, z1, x2, y2, z2)

#define CURVE vec_curve
