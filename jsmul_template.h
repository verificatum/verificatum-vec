
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

#include "jsmul_h_template.h"

void
FUNCTION_NAME(vec_jsmul_init, POSTFIX)
     (FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX) table,
      CURVE *curve,
      size_t len,
      size_t block_width)
{
  size_t i;             /* Index parameters. */
  size_t tab_len;       /* Size of a subtable. */

  /* We need the curve parameter in other implementations of
   * multiplication, so to use function pointers we need to keep it here
   * as well.
   */
  VEC_UNUSED(curve);

  /* Initialize length parameters. */
  table->len = len;
  table->block_width = block_width;
  if (len < block_width) {
    table->block_width = len;
  }
  table->tabs_len = (len + block_width - 1) / block_width;

  /* Allocate and initialize space for pointers to tables. */
  table->tabsx =
    (FIELD_ELEMENT_VAR **)malloc(table->tabs_len * sizeof(FIELD_ELEMENT_VAR *));
  table->tabsy =
    (FIELD_ELEMENT_VAR **)malloc(table->tabs_len * sizeof(FIELD_ELEMENT_VAR *));
  table->tabsz =
    (FIELD_ELEMENT_VAR **)malloc(table->tabs_len * sizeof(FIELD_ELEMENT_VAR *));

  tab_len = 1 << block_width;
  for (i = 0; i < table->tabs_len; i++)
    {

      /* Last block may be more narrow than the other. */
      if (i == table->tabs_len - 1
          && len - (table->tabs_len - 1) * block_width < block_width)
        {
          block_width = len - (table->tabs_len - 1) * block_width;
          tab_len = 1 << block_width;
        }

      /* Allocate and initialize a table. */
      table->tabsx[i] = ARRAY_MALLOC_INIT(tab_len);
      table->tabsy[i] = ARRAY_MALLOC_INIT(tab_len);
      table->tabsz[i] = ARRAY_MALLOC_INIT(tab_len);
    }
}

int
FUNCTION_NAME(vec_jsmul_clear, POSTFIX)
     (FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX) table)
{
  size_t i;                                /* Index variables. */
  size_t tabs_len = table->tabs_len;       /* Number of sub tables. */
  size_t block_width = table->block_width; /* Width of each sub table. */
  size_t tab_len;                          /* Size of each sub table. */

  tab_len = 1 << block_width;
  for (i = 0; i < tabs_len; i++)
    {

      /* Last block may have smaller width, but it is never zero. */
      if (i == tabs_len - 1)
        {
          block_width = table->len - (tabs_len - 1) * block_width;
          tab_len = 1 << block_width;
        }

      /* Deallocate table. */
      ARRAY_CLEAR_FREE(table->tabsx[i], tab_len);
      ARRAY_CLEAR_FREE(table->tabsy[i], tab_len);
      ARRAY_CLEAR_FREE(table->tabsz[i], tab_len);
    }

  /* Deallocate table of tables. */
  free(table->tabsx);
  free(table->tabsy);
  free(table->tabsz);

  /* Hack to avoid unused-variable warnings for the case where
     ARRAY_CLEAR_FREE expands to an expression not involving
     tab_len. */
  return tab_len;
}

void
FUNCTION_NAME(vec_jsmul_precomp, POSTFIX)
     (FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX) table,
      CURVE *curve,
      FIELD_ELEMENT_VAR *basesx, FIELD_ELEMENT_VAR *basesy,
      FIELD_ELEMENT_VAR *basesz)
{
  size_t i, j;                    /* Index variables. */
  size_t block_width;             /* Width of current subtable. */
  size_t tab_len;                 /* Size of current subtable. */
  FIELD_ELEMENT_VAR *tx;          /* Temporary variable for subtable. */
  FIELD_ELEMENT_VAR *ty;          /* Temporary variable for subtable. */
  FIELD_ELEMENT_VAR *tz;          /* Temporary variable for subtable. */

  int mask;           /* Mask used for dynamic programming */
  int one_mask;       /* Mask containing a single non-zero bit. */

  SCRATCH(scratch);

  VEC_UNUSED(curve);

  SCRATCH_INIT(scratch);

  block_width = table->block_width;
  tab_len = 1 << block_width;

  for (i = 0; i < table->tabs_len; i++)
    {

      /* Last block may have smaller width, but it is never zero. */
      if (i == table->tabs_len - 1)
        {
          block_width = table->len - (table->tabs_len - 1) * block_width;
          tab_len = 1 << block_width;
        }

      /* Current subtable. */
      tx = table->tabsx[i];
      ty = table->tabsy[i];
      tz = table->tabsz[i];

      /* Initialize current subtable with all trivial products. */
      FIELD_ELEMENT_VAR_UNIT(tx[0], ty[0], tz[0]);

      mask = 1;
      for (j = 0; j < block_width; j++)
        {
          FIELD_ELEMENT_VAR_SET(tx[mask], ty[mask], tz[mask],
                                basesx[j], basesy[j], basesz[j]);
          mask <<= 1;
        }

      /* Initialize current subtable with all non-trivial products. */
      for (mask = 1; ((size_t) mask) < tab_len; mask++)
        {
          one_mask = mask & (-mask);

          JADD_VAR(scratch,
                   tx[mask], ty[mask], tz[mask],
                   curve,
                   tx[mask ^ one_mask],
                   ty[mask ^ one_mask],
                   tz[mask ^ one_mask],
                   tx[one_mask],
                   ty[one_mask],
                   tz[one_mask]);
        }

      basesx += block_width;
      basesy += block_width;
      basesz += block_width;
    }

  SCRATCH_CLEAR(scratch);
}

static int
getbits(mpz_t *op, int index, size_t block_width)
{
  int i;
  int bits = 0;

  for (i = block_width - 1; i >= 0; i--)
    {
      bits <<= 1;
      if (mpz_tstbit(op[i], index))
        {
          bits |= 1;
        }
    }
  return bits;
}

void
FUNCTION_NAME(vec_jsmul_table, POSTFIX)
     (FIELD_ELEMENT_VAR ropx, FIELD_ELEMENT_VAR ropy, FIELD_ELEMENT_VAR ropz,
      CURVE *curve,
      FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX) table,
      mpz_t *scalars,
      size_t max_scalar_bitlen)
{

  size_t i;
  int index;
  int mask;

  mpz_t *scls;

  FIELD_ELEMENT tmpx;
  FIELD_ELEMENT tmpy;
  FIELD_ELEMENT tmpz;

  size_t len = table->len;
  size_t tabs_len = table->tabs_len;
  size_t block_width = table->block_width;
  size_t last_block_width = len - (tabs_len - 1) * block_width;

  FIELD_ELEMENT_VAR **tabsx = table->tabsx;
  FIELD_ELEMENT_VAR **tabsy = table->tabsy;
  FIELD_ELEMENT_VAR **tabsz = table->tabsz;

  SCRATCH(scratch);

  VEC_UNUSED(curve);

  SCRATCH_INIT(scratch);

  FIELD_ELEMENT_INIT(tmpx);
  FIELD_ELEMENT_INIT(tmpy);
  FIELD_ELEMENT_INIT(tmpz);

  /* Initialize result variable. */
  FIELD_ELEMENT_UNIT(tmpx, tmpy, tmpz);

  /* Execute simultaneous square-and-multiply. */
  for (index = max_scalar_bitlen - 1; index >= 0; index--)
    {

      /* Square ... */
      JDBL(scratch,
           tmpx, tmpy, tmpz,
           curve,
           tmpx, tmpy, tmpz);

      /* ... and multiply. */
      i = 0;
      scls = scalars;
      while (i < tabs_len)
        {
          if (i == tabs_len - 1)
            {
              mask = getbits(scls, index, last_block_width);
            }
          else
            {
              mask = getbits(scls, index, block_width);
            }

          JADD(scratch,
               tmpx, tmpy, tmpz,
               curve,
               tmpx, tmpy, tmpz,
               tabsx[i][mask], tabsy[i][mask], tabsz[i][mask]);

          i++;
          scls += block_width;
        }
    }

  SCRATCH_CLEAR(scratch);

  FIELD_ELEMENT_CONTRACT(ropx, ropy, ropz, tmpx, tmpy, tmpz);
  FIELD_ELEMENT_CLEAR(tmpx);
  FIELD_ELEMENT_CLEAR(tmpy);
  FIELD_ELEMENT_CLEAR(tmpz);
}

void
FUNCTION_NAME(vec_jsmul_block_batch, POSTFIX)
     (FIELD_ELEMENT ropx, FIELD_ELEMENT ropy, FIELD_ELEMENT ropz,
      CURVE *curve,
      FIELD_ELEMENT_VAR *basesx, FIELD_ELEMENT_VAR *basesy,
      FIELD_ELEMENT_VAR *basesz,
      mpz_t *scalars,
      size_t len,
      size_t block_width, size_t batch_len,
      size_t max_scalar_bitlen)
{
  size_t i;
  FUNCTION_NAME(vec_jsmul_tab, TAB_POSTFIX) table;

  FIELD_ELEMENT_VAR tmpx;
  FIELD_ELEMENT_VAR tmpy;
  FIELD_ELEMENT_VAR tmpz;

  SCRATCH(scratch);

  SCRATCH_INIT(scratch);

  FIELD_ELEMENT_VAR_INIT(tmpx);
  FIELD_ELEMENT_VAR_INIT(tmpy);
  FIELD_ELEMENT_VAR_INIT(tmpz);

  if (len < batch_len) {
    batch_len = len;
  }

  FUNCTION_NAME(vec_jsmul_init, POSTFIX)(table, curve, batch_len,
                                         block_width);

  /* Initialize result to unit element. */
  FIELD_ELEMENT_UNIT(ropx, ropy, ropz);

  for (i = 0; i < len; i += batch_len)
    {

      /* Last batch may be slightly shorter. */
      if (len - i < batch_len)
        {
          batch_len = len - i;

          FUNCTION_NAME(vec_jsmul_clear, POSTFIX)(table);
          FUNCTION_NAME(vec_jsmul_init, POSTFIX)(table, curve,
                                                 batch_len, block_width);
        }

      /* Perform computation for batch */
      FUNCTION_NAME(vec_jsmul_precomp, POSTFIX)(table, curve,
                                                basesx, basesy, basesz);

      /* Compute batch. */
      FUNCTION_NAME(vec_jsmul_table, POSTFIX)(tmpx, tmpy, tmpz,
                                              curve, table,
                                              scalars, max_scalar_bitlen);

      /* Add with result so far. */
      JADD(scratch,
           ropx, ropy, ropz,
           curve,
           ropx, ropy, ropz,
           tmpx, tmpy, tmpz);

      /* Move on to next batch. */
      basesx += batch_len;
      basesy += batch_len;
      basesz += batch_len;
      scalars += batch_len;
    }

  SCRATCH_CLEAR(scratch);

  FUNCTION_NAME(vec_jsmul_clear, POSTFIX)(table);
}

void
FUNCTION_NAME(vec_jsmul, POSTFIX)
     (FIELD_ELEMENT ropx, FIELD_ELEMENT ropy, FIELD_ELEMENT ropz,
      vec_curve *curve,
      FIELD_ELEMENT_VAR *basesx, FIELD_ELEMENT_VAR *basesy,
      FIELD_ELEMENT_VAR *basesz,
      mpz_t *scalars,
      size_t len)
{
  size_t i;
  size_t bitlen;
  size_t max_scalar_bitlen;
  size_t batch_len = 100;      /* This is somewhat arbitrary, but it
                                  makes the amortized cost for squaring
                                  very small in comparison to the cost
                                  for multiplications. */
  size_t block_width;

  /* Compute the maximal bit length among the scalars. */
  max_scalar_bitlen = 0;
  for (i = 0; i < len; i++)
    {
      bitlen = mpz_sizeinbase(scalars[i], 2);
      if (bitlen > max_scalar_bitlen)
        {
          max_scalar_bitlen = bitlen;
        }
    }

  /* Determine a good block width. */
  block_width = vec_smul_block_width(max_scalar_bitlen, batch_len);

  FUNCTION_NAME(vec_jsmul_block_batch, POSTFIX)(ropx, ropy, ropz,
                                                curve,
                                                basesx, basesy, basesz,
                                                scalars,
                                                len,
                                                block_width, batch_len,
                                                max_scalar_bitlen);
}
