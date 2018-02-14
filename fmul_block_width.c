
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

/*
 * Computes a "theoretical" optimal block width for a given scalar
 * length.
 */
int
vec_fmul_block_width(int bit_length, int len) {

  int width = 2;
  double cost = 1.5 * bit_length;
  double oldCost;
  double t;
  double s;
  double m;

  do {

    oldCost = cost;

    /* Amortized cost for table. */
    t = ((double)(1 << width) - width + bit_length) / len;

    /* Cost for doubling. */
    s = ((double)bit_length) / width;

    /* Cost for adding. */
    m = ((double)bit_length) / width;

    cost = t + m + s;

    width++;

  } while (width < 17 && cost < oldCost);

  width--;

  return width;
}
