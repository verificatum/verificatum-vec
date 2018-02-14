
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

static void mpz_t_to_felem(felem rop, const mpz_t op)
{
  memset(rop, 0, sizeof(felem));
  mpz_export (rop,           /* We write directly into the felem. */
              NULL,          /* We do not care how many bytes are copied. */
              -1,            /* Least significant GMP-limb of op goes
                                into the first limb of rop. */
              sizeof(limb),  /* Size of each limb. */
              0,             /* To native endianness. */
              8,             /* Most significant 8 bits should be zero. */
              op);
}

static void felem_to_mpz_t(mpz_t rop, const felem op)
{
  if (felem_is_zero(op))
    {
      mpz_set_si(rop, 0);
    }
  else
    {
      mpz_import(rop,
                 4,              /* Number of limbs in felem. */
                 -1,             /* First limb of op goes into the least
                                    significant limb of rop. */
                 sizeof(limb),   /* Size of each limb. */
                 0,              /* From native endianness. */
                 8,              /* Most significant 8 bits should be
                                    ignored */
                 op);
    }
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
