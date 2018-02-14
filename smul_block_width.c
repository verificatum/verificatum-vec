
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

/*
 * Computes a "theoretical" optimal block width for a given scalar
 * length.
 */
int
vec_smul_block_width(int scalars_bitlen, int batch_len) {

  /* This computes the theoretical optimum. */
  int width = 1;
  double cost = 1.5 * scalars_bitlen;
  double dbl;
  double add_precomp;
  double add;
  double old_cost;
  int width_exp;

  do {

    old_cost = cost;

    width++;
    width_exp = 1 << width;

    /* Amortized cost for doublings. */
    dbl = ((double)scalars_bitlen) / batch_len;

    /* Amortized cost for precomputing for a block. */
    add_precomp = ((double)width_exp) / width;

    /* Amortized cost for adding, given precomputation. */
    add = (((double)(scalars_bitlen)) / width) * (1 - 1.0 / width_exp);

    /* Total amortized cost. */
    cost = dbl + add_precomp + add;

  } while (cost < old_cost);

  width--;

  return width;
}
