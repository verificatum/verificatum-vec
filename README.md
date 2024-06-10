# Verificatum Elliptic Curve Library (VEC)

## Overview

This library implements elliptic curves based on the GMP library, for
simultaneous or fixed base modular exponentiation. The formulas used
are taken from https://www.hyperelliptic.org.

The purpose of this library is twofold. Firstly, we want a
self-contained and relatively small library that we hope can converge
to a static state. Secondly, our particular application allow using
certain optimizations and we need low-level access to these.

Our code is **not intended to be secure against side-channel attacks**,
since we do not need it for the applications we have in mind. It is
the responsibility of the user to make sure that this is the case in
their application.

For a detailed account of such algorithms, a good source is [Handbook
of Applied Cryptography](http://www.cacr.math.uwaterloo.ca/hac),
*Menezes, Oorshot, and Vanstone*, which is available for free. We use
additional techniques only suitable for very large number of
exponentiations.

We have copied optimized code from the [OpenSSL
project](https://www.openssl.org) for a few standard curves. These are
P-224 (written by Emilia Käsper), P-256 and P-521 (written by David
Langley). This code is in turn heavily inspired by the implementation
of Curve25519 by Dan Bernstein. The optimized code is roughly a factor
of three faster than the code based cleanly on top of GMP, which
explains the difference in running time between curves of the same
size.

Torbjorn Granlund helped us with the benchmarking. Emilia Käsper took
the time to answer our questions about her code and interpreting the
benchmarks. Dan Bernstein gave advice on how to implement the default
curves and pointed to Emilias code.

Some of the algorithms are implemented using macros. This may be
thought of as poor man's C++ templates. The rationale behind this
approach is that it is a way to include specialized add/double
routines implemented by others in such a way that their datatypes need
not be modified. Instead all that is needed is a macro file that maps
actual datatypes to macros.

This is necessary, since any conversion to/from a canonical type,
e.g., based on GMP, is too costly when the curve operations are
aggressively optimized. We avoid C++ to not involve yet another
language.

The following assumes that you are using a release. Developers should
also read `README_DEV.md`.


## Building

Building has been tested with GMP 6.1.2. Then `LIBRARY_PATH` must
point to `libgmp.la` and `C_INCLUDE_PATH` must point to `gmp.h`. This
is usually the case automatically after installing GMP.

Then use

        ./configure
        make

to build the library.


## Installing

Use

        make install

to install the library `libvec.{la,a,so}`.


You may need to run `sudo /sbin/ldconfig` on some platforms which have
flawed implementations of the cache that stores locations of
libraries.

If you prefer to use the [Clang compiler](https://clang.llvm.org) in
place of GCC for the native code, then you may use `./configure
CC=clang` instead of the above to enable it.

**Caution: Please understand that although it seems that Clang works
as well as GCC, switching compiler is a large change for mature
software.**


## Usage

The software is supposed to be used by other applications, but you can
use

        make check

or

        vec test

to test the arithmetic of all implemented curves, and you can use

        make bench

or

        vec speed

to get some benchmarks. Please consult the code to see exactly what is
measured before drawing any conclusions. Both commands also accept
names of curves as input to restrict the execution to these, e.g.,

        vec test P-224 P-256

tests the arithmetic of the curves P-224 and P-256 and nothing else.


## API Documentation

You may use

        make api

to build also some documentation using Doxygen (this assumes you have
installed `doxygen`). The API is not installed anywhere. You can copy
it to any location.


## Reporting Bugs

Minor bugs should be reported in the repository system as issues or
bugs. Security critical bugs, vulnerabilities, etc should be reported
directly to
[https://www.verificatum.org](https://www.verificatum.org). We will
make best effort to disclose the information in a responsible way
before the finder gets proper credit.
