
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

#if VERIFICATUM

static void
fprint_felem(FILE *out, char *name, felem f)
{
  int i;
  fprintf(out, "%s = ", name);
  for (i = 8; i >= 0; i--)
    {
      fprintf(out, "%04x", (unsigned int)(((limb)f[i]) >> 48) & 0xFFFF);
      fprintf(out, "%04x", (unsigned int)(((limb)f[i]) >> 32) & 0xFFFF);
      fprintf(out, "%04x", (unsigned int)(((limb)f[i]) >> 16) & 0xFFFF);
      fprintf(out, "%04x", (unsigned int)(f[i] & 0xFFFF));
      if (i > 0)
        {
          fprintf(out, "_");
        }
    }
  fprintf(out, "\n");
}

#endif /* VERIFICATUM */

static void felem_to_mpz_t(mpz_t rop, const felem op)
{
  felem cop;

  if (felem_is_zero(op))
    {
      mpz_set_si(rop, 0);
    }
  else
    {
      felem_contract(cop, op);
      mpz_import(rop,
                 9,              /* Number of limbs in felem. */
                 -1,             /* First limb of op goes into the least
                                    significant limb of rop. */
                 sizeof(limb),   /* Size of each limb. */
                 0,              /* From native endianness. */
                 6,              /* Most significant 6 bits should be
                                    ignored */
                 cop);
    }
}

static void mpz_t_to_felem(felem rop, const mpz_t op)
{
  memset(rop, 0, sizeof(felem));
  mpz_export (rop,           /* We write directly into the felem. */
              NULL,          /* We do not care how many bytes are copied. */
              -1,            /* Least significant GMP-limb of op goes
                                into the first limb of rop. */
              sizeof(limb),  /* Size of each limb. */
              0,             /* To native endianness. */
              6,             /* Most significant 6 bits should be zero. */
              op);
}

static felem* mpz_t_s_to_felems(mpz_t *ops, size_t len)
{
  size_t i;
  felem *res = (felem*)malloc(len * sizeof(felem));

  for (i = 0; i < len; i++)
    {
      mpz_t_to_felem(res[i], ops[i]);
    }
  return res;
}

#if VERIFICATUM

static mpz_t* felems_to_mpz_t_s(felem *ops, size_t len)
{
  int i;
  mpz_t *res = vec_array_alloc_init(len);

  for (i = 0; i < len; i++)
    {
      felem_to_mpz_t(res[i], ops[i]);
    }
  return res;
}

#endif /* VERIFICATUM */
