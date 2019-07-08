
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

/*! \mainpage README
 *
 * \verbinclude README
 */

/** @file */
#ifndef VEC_H
#define VEC_H

#include <gmp.h>

/*
 * ********************* BASIC ROUTINES *************************
 */

/**
 * Computes the modular square root of the input. The behaviour is
 * undefined for quadratic non-residues. This is used to embed
 * arbitrary strings into group elements.
 */
void
vec_sqrt(mpz_t res, mpz_t a, mpz_t p);

/**
 * Simple alarm predicate.
 *
 * @param start_time Point in time the clock started ticking.
 * @param interval For how long the clock should be ticking.
 */
int
vec_done(long start_time, long interval);

/**
 * Allocates an array of uninitialized mpz_t instances.
 *
 * @param len Number of instances.
 */
mpz_t *
vec_array_alloc(size_t len);

/**
 * Allocates an array of initialized mpz_t instances.
 *
 * @param len Number of instances.
 */
mpz_t *
vec_array_alloc_init(size_t len);

/**
 * Clears all instances in the input array and then frees the array.
 *
 * @param a Array of mpz_t instances.
 * @param len Number of instances.
 */
void
vec_array_clear_free(mpz_t *a, size_t len);


/**
 * Block of temporary variables used by doubling and adding algorithms
 * for the curves.
 */
typedef struct {
  mpz_t t1; /**< Temorary scratch variable. */
  mpz_t t2; /**< Temorary scratch variable. */
  mpz_t t3; /**< Temorary scratch variable. */
  mpz_t t4; /**< Temorary scratch variable. */
  mpz_t t5; /**< Temorary scratch variable. */
  mpz_t t6; /**< Temorary scratch variable. */
  mpz_t t7; /**< Temorary scratch variable. */
  mpz_t t8; /**< Temorary scratch variable. */
  mpz_t t9; /**< Temorary scratch variable. */
} vec_scratch_mpz_t[1];

/**
 * Initalizes the mpz_t instances of the block.
 *
 * @param scratch Block of temporary variables.
 */
void
vec_scratch_init_mpz_t(vec_scratch_mpz_t scratch);

/**
 * Clears the mpz_t instances of the block.
 *
 * @param scratch Block of temporary variables.
 */
void
vec_scratch_clear_mpz_t(vec_scratch_mpz_t scratch);


/*
 * ************** TYPES FOR CURVE ALGORITHMS *********************
 */

/**
 * Type for timer function for the core doubling and adding
 * algorithms.
 */
typedef long (*coretimer_func)(int test_time, mpz_t X, mpz_t Y);


/**
 * Struct representing a curve including the algorithms used to
 * compute in the curve.
 */
struct vec_curve;

/**
 * Union "pointer" to distinct structs. This is convenient when
 * passing pointers over a Java Native Interface (JNI).
 */
typedef union
{
  struct _vec_jfmul_tab_generic_inner *generic;   /**< Generic table. */
  struct _vec_jfmul_tab_nistp224_inner *nistp224; /**< nistp224 table. */
  struct _vec_jfmul_tab_nistp256_inner *nistp256; /**< nistp256 table. */
  struct _vec_jfmul_tab_nistp521_inner *nistp521; /**< nistp521 table. */

} vec_jfmul_tab_ptr;

/**
 * Doubling algorithm using Jacobi coordinates.
 */
typedef void (*jdbl_func)(vec_scratch_mpz_t scratch,
                          mpz_t X3, mpz_t Y3, mpz_t Z3,
                          struct vec_curve *curve,
                          mpz_t X1, mpz_t Y1, mpz_t Z1);

/**
 * Addition algorithm using Jacobi coordinates.
 */
typedef void (*jadd_func)(vec_scratch_mpz_t scratch,
                          mpz_t X3, mpz_t Y3, mpz_t Z3,
                          struct vec_curve *curve,
                          mpz_t X1, mpz_t Y1, mpz_t Z1,
                          mpz_t X2, mpz_t Y2, mpz_t Z2);

/**
 * Multiplication algorithm using Jacobi coordinates.
 */
typedef void (*jmul_func)(mpz_t RX, mpz_t RY, mpz_t RZ,
                          struct vec_curve *curve,
                          mpz_t X, mpz_t Y, mpz_t Z,
                          mpz_t scalar);

/**
 * Simultaneous multiplication algorithm using Jacobi coordinates.
 */
typedef void (*jsmul_func)(mpz_t RX, mpz_t RY, mpz_t RZ,
                           struct vec_curve *curve,
                           mpz_t *X, mpz_t *Y, mpz_t *Z,
                           mpz_t *scalars,
                           size_t len);

/**
 * Precomputation for fixed basis multiplication algorithm using
 * Jacobi coordinates.
 */
typedef vec_jfmul_tab_ptr (*jfmul_precomp_func)(struct vec_curve *curve,
                                                mpz_t X, mpz_t Y, mpz_t Z,
                                                size_t len);

/**
 * Algorithm for fixed basis multiplication algorithm using Jacobi
 * coordinates.
 */
typedef void (*jfmul_func)(mpz_t RX, mpz_t RY, mpz_t RZ,
                           struct vec_curve *curve,
                           vec_jfmul_tab_ptr ptr,
                           mpz_t scalar);

/**
 * Algorithm for freeing resources allocated during precomputation for
 * fixed basis multiplication.
 */
typedef void (*jfmul_free_func)(vec_jfmul_tab_ptr table);


/*
 * ********************* CURVE MANIPULATION *************************
 */


/**
 * Struct representing a curve including the algorithms used to
 * compute in the curve.
 */
struct vec_curve {
  const char *name;                  /**< Standard name of curve. */
  mpz_t modulus;                     /**< Prime modulus. */
  mpz_t a;                           /**< x-coefficient. */
  mpz_t b;                           /**< Constant coefficient. */
  mpz_t gx;                          /**< x-coefficient of generator. */
  mpz_t gy;                          /**< x-coefficient of generator. */
  mpz_t n;                           /**< Order of curve. */

  jdbl_func jdbl;                    /**< Doubling function. */
  jadd_func jadd;                    /**< Addition function. */
  jmul_func jmul;                    /**< Multiplication function. */
  jsmul_func jsmul;                  /**< Simultaneous multiplication
                                        function. */
  jfmul_precomp_func jfmul_precomp;  /**< Fixed base pre-computation function.*/
  jfmul_func jfmul;                  /**< Fixed base multiplication function.*/
  jfmul_free_func jfmul_free;        /**< Free fixed base table function.*/
  coretimer_func jdbl_timer;         /**< Timer function for doubling.*/
  coretimer_func jadd_timer;         /**< Timer function for addition.*/
};

/**
 * Struct representing a curve including the algorithms used to
 * compute in the curve.
 */
typedef struct vec_curve vec_curve;

/**
 * Allocates a curve.
 */
vec_curve *
vec_curve_alloc();

/**
 * Frees a curve.
 */
void
vec_curve_free(vec_curve *curve);

/**
 * Returns a curve with the given explicit curve parameters. The
 * default algorithms are used.
 */
vec_curve *
vec_curve_get_anon(mpz_t modulus, mpz_t a, mpz_t b,
                   mpz_t gx, mpz_t gy, mpz_t n);

/**
 * Returns the named curve. Optimized algorithms are used if
 * available.
 */
vec_curve *
vec_curve_get_named(char *name, int implementation);

/**
 * Returns the named curve. Optimized algorithms are used if
 * available. Name is a string without a NULL at the end.
 */
vec_curve *
vec_curve_get_named_len(char *name, int len, int implementation);

/**
 * Returns the number of named curve. This is used to iterate through
 * all curves, e.g. for testing.
 */
int
vec_curve_number_of_names();

/**
 * Returns the name of the ith curve.
 */
char *
vec_curve_get_name(int i);

/**
 * Predicate for the property that a = -3, where a is the linear curve
 * coefficient.
 */
int
vec_curve_a_eq_neg3(vec_curve *curve);

/**
 * Equality predicate for points in affine coordinates.
 */
int
vec_eq(mpz_t x1, mpz_t y1, mpz_t x2, mpz_t y2);




/*******************************************************************
 ** ARITHMETIC FOR CURVES IN AFFINE COORDINATES USED FOR DEBUGGING *
 *******************************************************************/

/**
 * Computes the optimal block width to be used during simultaneous
 * multiplication.
 */
int
vec_smul_block_width(int scalars_bitlen, int batch_len);

/**
 * Returns the optimal block width for fixed basis multiplication.
 */
int
vec_fmul_block_width(int bit_length, int len);

/**
 * Computes the doubling of the input point in affine coordinates.
 */
void
vec_dbl(vec_scratch_mpz_t scratch,
        mpz_t rx, mpz_t ry,
        vec_curve *curve,
        mpz_t x, mpz_t y);

/**
 * Computes the sum of the two input points in affine coordinates.
 */
void
vec_add(vec_scratch_mpz_t scratch,
        mpz_t rx, mpz_t ry,
        vec_curve *curve,
        mpz_t x1, mpz_t y1,
        mpz_t x2, mpz_t y2);

/**
 * Compute the scalar multiple of the input point in affine
 * coordinates.
 */
void
vec_mul(mpz_t rx, mpz_t ry,
        vec_curve *curve,
        mpz_t x, mpz_t y,
        mpz_t scalar);

/**
 * Table storing precomputed values during simultaneous
 * multiplication.
 */
typedef struct
{
  vec_curve *curve;       /**< Underlying curve. */
  size_t len;             /**< Total number of bases/scalars. */
  size_t block_width;     /**< Number of bases/scalars in each block. */
  size_t tabs_len;        /**< Number of blocks. */
  mpz_t **tabsx;          /**< Table of tables, one sub-table for each block. */
  mpz_t **tabsy;          /**< Table of tables, one sub-table for each block. */

} vec_smul_tab[1]; /* Magic references. */


/**
 * Initializes the table for simultaneous multiplication.
 */
void
vec_smul_init(vec_smul_tab table,
              vec_curve *curve,
              size_t len, size_t block_width);

/**
 * Releases memory allocated for the table for simultaneous
 * multiplication.
 */
void
vec_smul_clear(vec_smul_tab table);

/**
 * Performs pre-computation for the given table and bases.
 */
void
vec_smul_precomp(vec_smul_tab table, mpz_t *basesx, mpz_t *basesy);

/**
 * Raises the bases of the pre-computed table simultaneously to the
 * given scalars.
 */
void
vec_smul_table(mpz_t ropx, mpz_t ropy,
               vec_smul_tab table,
               mpz_t *scalars,
               size_t max_scalar_bitlen);

/**
 * Raises the bases of the pre-computed table simultaneously to the
 * given mulonents. This is done in batches of the given size.
 */
void
vec_smul_block_batch(mpz_t ropx, mpz_t ropy,
                     vec_curve *curve,
                     mpz_t *basesx, mpz_t *basesy,
                     mpz_t *scalars,
                     size_t len,
                     size_t block_width, size_t batch_len,
                     size_t max_scalar_bitlen);

/**
 * Computes the simultaneous multiplication of the points and scalars.
 */
void
vec_smul(mpz_t ropx, mpz_t ropy,
         vec_curve *curve,
         mpz_t *basesx, mpz_t *basesy,
         mpz_t *scalars,
         size_t len);



/*******************************************************************
 ******** ARITHMETIC FOR GENERAL CURVES IN JACOBI COORDINATES ******
 *******************************************************************/

/**
 * Computes the doubling of the input point in Jacobi coordinates.
 */
void
vec_jdbl_generic(vec_scratch_mpz_t scratch,
                 mpz_t X3, mpz_t Y3, mpz_t Z3,
                 vec_curve *curve,
                 mpz_t X1, mpz_t Y1, mpz_t Z1);

/**
 * Computes the sum of the two input points in Jacobi coordinates.
 */
void
vec_jadd_generic(vec_scratch_mpz_t scratch,
                 mpz_t X3, mpz_t Y3, mpz_t Z3,
                 vec_curve *curve,
                 mpz_t X1, mpz_t Y1, mpz_t Z1,
                 mpz_t X2, mpz_t Y2, mpz_t Z2);

/* Naive version of multiplication. Only used during development.
void
vec_jmul_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                 vec_curve *curve,
                 mpz_t X, mpz_t Y, mpz_t Z,
                 mpz_t scalar);
*/

/**
 * Compute the scalar multiple of the input point in Jacobi
 * coordinates using sliding window.
 */
void
vec_jmulsw_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                   vec_curve *curve,
                   mpz_t X, mpz_t Y, mpz_t Z,
                   mpz_t scalar);


/**
 * Computes the simultaneous multiplication of the points and scalars.
 */
void
vec_jsmul_generic(mpz_t ropx, mpz_t ropy, mpz_t ropz,
                  vec_curve *curve,
                  mpz_t *basesx, mpz_t *basesy, mpz_t *basesz,
                  mpz_t *scalars,
                  size_t len);

/**
 * Computes the doubling of the input point in Jacobi coordinates.
 */
void
vec_jdbl_a_eq_neg3_generic(vec_scratch_mpz_t scratch,
                           mpz_t X3, mpz_t Y3, mpz_t Z3,
                           vec_curve *curve,
                           mpz_t X1, mpz_t Y1, mpz_t Z1);

/**
 * Computes the sum of the two input points in Jacobi coordinates.
 */
void
vec_jadd_a_eq_neg3_generic(vec_scratch_mpz_t scratch,
                           mpz_t X3, mpz_t Y3, mpz_t Z3,
                           vec_curve *curve,
                           mpz_t X1, mpz_t Y1, mpz_t Z1,
                           mpz_t X2, mpz_t Y2, mpz_t Z2);

/* Naive version of multiplication. Only used during development.
void
vec_jmul_a_eq_neg3_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                           vec_curve *curve,
                           mpz_t X, mpz_t Y, mpz_t Z,
                           mpz_t scalar);
*/

/**
 * Compute the scalar multiple of the input point in Jacobi
 * coordinates using sliding window.
 */
void
vec_jmulsw_a_eq_neg3_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                             vec_curve *curve,
                             mpz_t X, mpz_t Y, mpz_t Z,
                             mpz_t scalar);

/**
 * Computes the simultaneous multiplication of the points and scalars.
 */
void
vec_jsmul_a_eq_neg3_generic(mpz_t ropx, mpz_t ropy, mpz_t ropz,
                            vec_curve *curve,
                            mpz_t *basesx, mpz_t *basesy, mpz_t *basesz,
                            mpz_t *scalars,
                            size_t len);




/*******************************************************************
 ***** ARITHMETIC FOR CURVES WITH a = -3 IN JACOBI COORDINATES *****
 *******************************************************************/

/**
 * Computes the doubling of the input point in Jacobi coordinates.
 */
void
vec_jdbl_a_eq_neg3_generic(vec_scratch_mpz_t scratch,
                   mpz_t X3, mpz_t Y3, mpz_t Z3,
                   vec_curve *curve,
                   mpz_t X1, mpz_t Y1, mpz_t Z1);

/**
 * Compute the scalar multiple of the input point in Jacobi
 * coordinates using double-and-add.
 */
void
vec_jmul_a_eq_neg3_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                   vec_curve *curve,
                   mpz_t X, mpz_t Y, mpz_t Z,
                   mpz_t scalar);

/**
 * Compute the scalar multiple of the input point in Jacobi
 * coordinates using sliding window.
 */
void
vec_jmulsw_a_eq_neg3_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                     vec_curve *curve,
                     mpz_t X, mpz_t Y, mpz_t Z,
                     mpz_t scalar);
/**
 * Computes the simultaneous multiplication of the points and scalars.
 */
void
vec_jsmul_a_eq_neg3_generic(mpz_t ropx, mpz_t ropy, mpz_t ropz,
                  vec_curve *curve,
                  mpz_t *basesx, mpz_t *basesy, mpz_t *basesz,
                  mpz_t *scalars,
                  size_t len);

/**
 * Performs precomputation for fixed basis multiplication in Jacobi
 * coordinates.
 */
vec_jfmul_tab_ptr
vec_jfmul_precomp_generic(vec_curve *curve,
                          mpz_t X, mpz_t Y, mpz_t Z,
                          size_t len);

/**
 * Computes a fixed basis multiplication in Jacobi coordinates.
 */
void
vec_jfmul_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                  vec_curve *curve, vec_jfmul_tab_ptr table,
                  mpz_t scalar);

/**
 * Frees allocated memory for fixed basis multiplication in Jacobi
 * coordinates.
 */
void
vec_jfmul_free_generic(vec_jfmul_tab_ptr ptr);

/**
 * Performs precomputation for fixed basis multiplication in Jacobi
 * coordinates.
 */
vec_jfmul_tab_ptr
vec_jfmul_precomp_a_eq_neg3_generic(vec_curve *curve,
                                    mpz_t X, mpz_t Y, mpz_t Z,
                                    size_t len);

/**
 * Computes a fixed basis multiplication in Jacobi coordinates.
 */
void
vec_jfmul_a_eq_neg3_generic(mpz_t RX, mpz_t RY, mpz_t RZ,
                            vec_curve *curve, vec_jfmul_tab_ptr table,
                            mpz_t scalar);

/**
 * Frees allocated memory for fixed basis multiplication in Jacobi
 * coordinates.
 */
void
vec_jfmul_free_a_eq_neg3_generic(vec_jfmul_tab_ptr ptr);

/**
 * Transforms the point to the standard affine form, i.e., Z=1 and X
 * and Y positive, or X and Y are both -1 to indicate the point at
 * infinity.
 */
void
vec_jaff(mpz_t X, mpz_t Y, mpz_t Z, vec_curve *curve);

/**
 * Transforms the point to the standard Jacobi form, i.e., Z=1 or Z=0
 * (in the case of the unit point input).
 */
void
vec_affj(mpz_t X, mpz_t Y, mpz_t Z);

/**
 * Doubles the input point using Jacobi coordinates internally and
 * then converts the result to affine coordinates.
 */
void
vec_jdbl_aff(vec_scratch_mpz_t scratch,
             mpz_t rx, mpz_t ry,
             vec_curve *curve,
             mpz_t x, mpz_t y);

/**
 * Computes the sum of the input points using Jacobi coordinates
 * internally and then converts the result to affine coordinates.
 */
void
vec_jadd_aff(vec_scratch_mpz_t scratch,
             mpz_t rx, mpz_t ry,
             vec_curve *curve,
             mpz_t x1, mpz_t y1,
             mpz_t x2, mpz_t y2);

/**
 * Compute the scalar multiple of the input point using Jacobi
 * coordinates internally and then converts the result to affine
 * coordinates.
 */
void
vec_jmul_aff(mpz_t rx, mpz_t ry,
             vec_curve *curve,
             mpz_t x, mpz_t y,
             mpz_t scalar);

/**
 * Compute the simultaneous multiplication of the input points and
 * scalars using Jacobi coordinates internally and then converts the
 * result to affine coordinates.
 */
void
vec_jsmul_aff(mpz_t ropx, mpz_t ropy,
              vec_curve *curve,
              mpz_t *basesx, mpz_t *basesy,
              mpz_t *exponents,
              size_t len);

/**
 * Perform precomputation for fixed basis multiplication using Jacobi
 * coordinates internally.
 */
vec_jfmul_tab_ptr
vec_jfmul_precomp_aff(vec_curve *curve,
                      mpz_t x, mpz_t y,
                      size_t len);

/**
 * Compute the fixed basis scalar multiple using Jacobi coordinates
 * internally and then converts the result to affine coordinates.
 */
void
vec_jfmul_aff(mpz_t rx, mpz_t ry,
              vec_curve *curve,
              vec_jfmul_tab_ptr table,
              mpz_t scalar);

/**
 * Frees the memory allocated for the table for fixed basis
 * multiplication.
 */
void
vec_jfmul_free_aff(vec_curve *curve, vec_jfmul_tab_ptr ptr);




/*******************************************************************
 ***** OPTIMIZED ARITHMETIC FOR CURVES IN JACOBI COORDINATES *******
 *******************************************************************/

/*
 * Emilia Käsper's implementation of nistp256/P-256.
 */

/*! @copydoc vec_jdbl_generic() */
void
vec_jdbl_nistp224(vec_scratch_mpz_t scratch,
                  mpz_t X3, mpz_t Y3, mpz_t Z3,
                  vec_curve *curve,
                  mpz_t X1, mpz_t Y1, mpz_t Z1);

/*! @copydoc vec_jadd_generic() */
void
vec_jadd_nistp224(vec_scratch_mpz_t scratch,
                  mpz_t X3, mpz_t Y3, mpz_t Z3,
                  vec_curve *curve,
                  mpz_t X1, mpz_t Y1, mpz_t Z1,
                  mpz_t X2, mpz_t Y2, mpz_t Z2);

/* Naive version of multiplication. Only used during development.
void
vec_jmul_nistp224(mpz_t RX, mpz_t RY, mpz_t RZ,
                  vec_curve *curve,
                  mpz_t X, mpz_t Y, mpz_t Z,
                  mpz_t scalar);
*/

/*! @copydoc vec_jmulsw_generic() */
void
vec_jmulsw_nistp224(mpz_t RX, mpz_t RY, mpz_t RZ,
                    vec_curve *curve,
                    mpz_t X, mpz_t Y, mpz_t Z,
                    mpz_t scalar);

/*! @copydoc vec_jsmul_generic() */
void
vec_jsmul_nistp224(mpz_t ropx, mpz_t ropy, mpz_t ropz,
                   vec_curve *curve,
                   mpz_t *basesx, mpz_t *basesy, mpz_t *basesz,
                   mpz_t *scalars,
                   size_t len);

/*! @copydoc vec_jfmul_precomp_generic() */
vec_jfmul_tab_ptr
vec_jfmul_precomp_nistp224(vec_curve *curve,
                           mpz_t X, mpz_t Y, mpz_t Z,
                           size_t len);

/*! @copydoc vec_jfmul_generic() */
void
vec_jfmul_nistp224(mpz_t RX, mpz_t RY, mpz_t RZ,
                   vec_curve *curve, vec_jfmul_tab_ptr ptr,
                   mpz_t scalar);

/*! @copydoc vec_jfmul_free_generic() */
void
vec_jfmul_free_nistp224(vec_jfmul_tab_ptr ptr);



/*
 * Adam Langley's implementation of nistp256/P-256.
 */

/*! @copydoc vec_jdbl_generic() */
void
vec_jdbl_nistp256(vec_scratch_mpz_t scratch,
                  mpz_t X3, mpz_t Y3, mpz_t Z3,
                  vec_curve *curve,
                  mpz_t X1, mpz_t Y1, mpz_t Z1);

/*! @copydoc vec_jadd_generic() */
void
vec_jadd_nistp256(vec_scratch_mpz_t scratch,
                  mpz_t X3, mpz_t Y3, mpz_t Z3,
                  vec_curve *curve,
                  mpz_t X1, mpz_t Y1, mpz_t Z1,
                  mpz_t X2, mpz_t Y2, mpz_t Z2);

/* Naive version of multiplication. Only used during development.
void
vec_jmul_nistp256(mpz_t RX, mpz_t RY, mpz_t RZ,
                  vec_curve *curve,
                  mpz_t X, mpz_t Y, mpz_t Z,
                  mpz_t scalar);
*/

/*! @copydoc vec_jmulsw_generic() */
void
vec_jmulsw_nistp256(mpz_t RX, mpz_t RY, mpz_t RZ,
                    vec_curve *curve,
                    mpz_t X, mpz_t Y, mpz_t Z,
                    mpz_t scalar);

/*! @copydoc vec_jsmul_generic() */
void
vec_jsmul_nistp256(mpz_t ropx, mpz_t ropy, mpz_t ropz,
                   vec_curve *curve,
                   mpz_t *basesx, mpz_t *basesy, mpz_t *basesz,
                   mpz_t *scalars,
                   size_t len);

/*! @copydoc vec_jfmul_precomp_generic() */
vec_jfmul_tab_ptr
vec_jfmul_precomp_nistp256(vec_curve *curve,
                           mpz_t X, mpz_t Y, mpz_t Z,
                           size_t len);

/*! @copydoc vec_jfmul_generic() */
void
vec_jfmul_nistp256(mpz_t RX, mpz_t RY, mpz_t RZ,
                   vec_curve *curve, vec_jfmul_tab_ptr ptr,
                   mpz_t scalar);

/*! @copydoc vec_jfmul_free_generic() */
void
vec_jfmul_free_nistp256(vec_jfmul_tab_ptr ptr);


/*
 * Adam Langley's implementation of nistp521/P-521.
 */

/*! @copydoc vec_jdbl_generic() */
void
vec_jdbl_nistp521(vec_scratch_mpz_t scratch,
                  mpz_t X3, mpz_t Y3, mpz_t Z3,
                  vec_curve *curve,
                  mpz_t X1, mpz_t Y1, mpz_t Z1);

/*! @copydoc vec_jadd_generic() */
void
vec_jadd_nistp521(vec_scratch_mpz_t scratch,
                  mpz_t X3, mpz_t Y3, mpz_t Z3,
                  vec_curve *curve,
                  mpz_t X1, mpz_t Y1, mpz_t Z1,
                  mpz_t X2, mpz_t Y2, mpz_t Z2);

/* Naive version of multiplication. Only used during development.
void
vec_jmul_nistp521(mpz_t RX, mpz_t RY, mpz_t RZ,
                  vec_curve *curve,
                  mpz_t X, mpz_t Y, mpz_t Z,
                  mpz_t scalar);
*/

/*! @copydoc vec_jmulsw_generic() */
void
vec_jmulsw_nistp521(mpz_t RX, mpz_t RY, mpz_t RZ,
                    vec_curve *curve,
                    mpz_t X, mpz_t Y, mpz_t Z,
                    mpz_t scalar);

/*! @copydoc vec_jsmul_generic() */
void
vec_jsmul_nistp521(mpz_t ropx, mpz_t ropy, mpz_t ropz,
                   vec_curve *curve,
                   mpz_t *basesx, mpz_t *basesy, mpz_t *basesz,
                   mpz_t *scalars,
                   size_t len);

/*! @copydoc vec_jfmul_precomp_generic() */
vec_jfmul_tab_ptr
vec_jfmul_precomp_nistp521(vec_curve *curve,
                           mpz_t X, mpz_t Y, mpz_t Z,
                           size_t len);

/*! @copydoc vec_jfmul_generic() */
void
vec_jfmul_nistp521(mpz_t RX, mpz_t RY, mpz_t RZ,
                   vec_curve *curve, vec_jfmul_tab_ptr ptr,
                   mpz_t scalar);

/*! @copydoc vec_jfmul_free_generic() */
void
vec_jfmul_free_nistp521(vec_jfmul_tab_ptr ptr);




/*
 * **** TIMING FUNCTIONS ********
 */

/**
 * Time doubling in Jacobi coordinates.
 */
long
time_jdbl_generic(int test_time, mpz_t X, mpz_t Y);

/**
 * Time addition in Jacobi coordinates.
 */
long
time_jadd_generic(int test_time, mpz_t X, mpz_t Y);


/**
 * Time doubling in Jacobi coordinates for P-224.
 */
long
time_jdbl_nistp224(int test_time, mpz_t X, mpz_t Y);

/**
 * Time addition in Jacobi coordinates for P-224.
 */
long
time_jadd_nistp224(int test_time, mpz_t X, mpz_t Y);


/**
 * Time doubling in Jacobi coordinates for P-256.
 */
long
time_jdbl_nistp256(int test_time, mpz_t X, mpz_t Y);

/**
 * Time addition in Jacobi coordinates for P-256.
 */
long
time_jadd_nistp256(int test_time, mpz_t X, mpz_t Y);


/**
 * Time doubling in Jacobi coordinates for P-521.
 */
long
time_jdbl_nistp521(int test_time, mpz_t X, mpz_t Y);

/**
 * Time addition in Jacobi coordinates for P-521.
 */
long
time_jadd_nistp521(int test_time, mpz_t X, mpz_t Y);


#endif /* VEC_H */
