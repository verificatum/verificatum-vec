
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

#include <stdio.h>
#include <gmp.h>
#include "vec.h"

#include "generic_macros.h"
#include "jfmul_template.h"

size_t
vec_jfmul_ptrsize()
{
  return sizeof(vec_jfmul_tab*);
}

vec_jfmul_tab_ptr
vec_jfmul_precomp(vec_curve *curve,
                    mpz_t X, mpz_t Y, mpz_t Z,
                    size_t len)
{
  vec_jfmul_tab_ptr ptr;

  ptr.generic = (vec_jfmul_tab*)malloc(sizeof(vec_jfmul_tab));

  vec_jfmul_init(ptr.generic, curve, len);

  vec_jfmul_prcmp(curve, ptr.generic, X, Y, Z);

  return ptr;
}

void
vec_jfmul(mpz_t RX, mpz_t RY, mpz_t RZ,
            vec_curve *curve, vec_jfmul_tab_ptr table_ptr,
            mpz_t scalar)
{

  vec_jfmul_cmp(RX, RY, RZ,
                  curve, table_ptr.generic,
                  scalar);
}

void
vec_jfmul_free(vec_jfmul_tab_ptr table_ptr)
{
  vec_jfmul_clear_free(table_ptr.generic);
}
