
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

#include "generic_macros.h"

#undef POSTFIX

#define POSTFIX _a_eq_neg3_generic_inner
#define TAB_POSTFIX _generic_inner

#undef JDBL

#define JDBL(scratch, rx, ry, rz, curve, x, y, z) \
  vec_jdbl_a_eq_neg3_generic(scratch,             \
                             rx, ry, rz,          \
                             curve,               \
                             x, y, z)
