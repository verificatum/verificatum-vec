
# Copyright 2008-2018 Douglas Wikstrom
#
# This file is part of Verificatum Elliptic Curve library (VEC).
#
# VEC is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# VEC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General
# Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with VEC. If not, see <http://www.gnu.org/licenses/>.

.PHONY: clean cleanapi

ACLOCAL_AMFLAGS=-I m4

BINDIR = bin

# Figuring out the best compiler/building flags to use for a given
# platform is difficult. This project uses GMP to implement large
# integers, so the flags used by GMP to compile should be a good
# default on all our target platforms.
#
# The authors of GMP had the good taste to provide these flags as the
# macro __GMP_CFLAGS, so writing a program to print the flags is
# trivial.

GMP_CFLAGS := $(shell ./extract_GMP_CFLAGS)

OPTIONAL_FLAGS=-O3
# OPTIONAL_FLAGS=-fprofile-arcs -ftest-coverage

# We use pedantic flags, strip the optimization flag of GMP, and
# insert our own level of optimization.
AM_CFLAGS := -Wall -W -Werror $(shell echo ${GMP_CFLAGS} | sed -e "s/-O[O12345]//") $(OPTIONAL_FLAGS)

AM_LDFLAGS = -lgmp

scriptmacros.m4:
	@printf "define(M4_VERSION, $(VERSION))dnl\n" > scriptmacros.m4

$(BINDIR)/vec-info: scriptmacros.m4 vec-info.src
	mkdir -p $(BINDIR)
	cat scriptmacros.m4 vec-info.src | m4 > $(BINDIR)/vec-info
	chmod +x $(BINDIR)/vec-info

UTILITY_SOURCES = array_alloc.c array_alloc_init.c array_clear_free.c done.c
MPZ_T_SOURCES = scratch_init_mpz_t.c scratch_clear_mpz_t.c
TABLE_OPTIMIZE_SOURCES = smul_block_width.c fmul_block_width.c
NAIVE_SOURCES = dbl.c add.c mul.c smul_init.c smul_clear.c smul_precomp.c smul_table.c smul_block_batch.c smul.c
GENERIC_SOURCES = jdbl_generic_inner.c jdbl_a_eq_neg3_generic_inner.c jadd_generic_inner.c
INNER_SOURCES = generic.c a_eq_neg3_generic.c nistp224.c nistp256.c nistp521.c 
AFFINE_SOURCES = jfmul_precomp_aff.c jfmul_aff.c jfmul_free_aff.c jaff.c affj.c jdbl_aff.c jadd_aff.c jmul_aff.c jsmul_aff.c
CURVE_SOURCES = curve_alloc.c curve_free.c curve_get_named.c eq.c sqrt.c

lib_LTLIBRARIES = libvec.la
libvec_la_SOURCES = ${UTILITY_SOURCES} ${MPZ_T_SOURCES} ${TABLE_OPTIMIZE_SOURCES} ${NAIVE_SOURCES} ${GENERIC_SOURCES} ${INNER_SOURCES} ${AFFINE_SOURCES} ${CURVE_SOURCES}

libvec_la_LIBADD = -lgmp
vec_LDADD = -lvec

include_HEADERS = vec.h
vec_SOURCES = vec.c vec.h
bin_PROGRAMS = vec
dist_bin = $(BINDIR)/vec-info
dist_bin_SCRIPTS = $(BINDIR)/vec-info

dist_noinst_DATA = extract_GMP_CFLAGS.c README.md COPYING NEWS AUTHORS ChangeLog config.h jmul_template.h nistp224_macros.h vec.h jsmul_h_template.h nistp256_macros.h jfmul_h_template.h jsmul_template.h nistp521_macros.h jfmul_template.h templates.h jmulsw_template.h generic_macros.h a_eq_neg3_generic_macros.h undefine_macros.h ecp_nistp224_core.c ecp_nistp256_core.c ecp_nistp521_core.c ecp_nistp224_util.c ecp_nistp256_util.c ecp_nistp521_util.c doxygen.cfg vec-info.src

all-local: .check_info.stamp

.check_info.stamp:
	@echo ""
	@echo "================================================================"
	@echo ""
	@echo "          INSTALLING AND VERIFYING INSTALLATION"
	@echo ""
	@echo " To install, you now use the standard 'make install' (or rather"
	@echo " 'sudo make install' on most Linuxes), but before you do that "
	@echo " it is prudent to use 'make check' to verify the binary."
	@echo " Please be patient, since this verifies the arithmetic of each"
	@echo " curve separately."
	@echo ""
	@echo "================================================================"
	@echo ""
	@touch .check_info.stamp

install-data-hook:
	@echo ""
	@echo "================================================================"
	@echo ""
	@echo "          TESTING AND BENCHMARKING"
	@echo ""
	@echo " Please consult the README file for information about how to"
	@echo " use the vec executable to test and measure the speed of the"
	@echo " library, or simply execute 'vec' to print usage information."
	@echo ""
	@echo "================================================================"
	@echo ""

check: vec
	@./vec check

coverage: all
	@./vec check
	@rm -rf coverage
	@mkdir -p coverage
	@lcov --capture --directory . --output-file coverage/coverage.info
	@genhtml coverage/coverage.info --output-directory coverage/html

bench: vec
	@./vec speed

api:
	@mkdir -p html_api
	doxygen doxygen.cfg 
	mv html_api/html api
	rm -rf html_api

cleancoverage:
	@rm -rf *.gcno *.gcda coverage coverage.info

cleanapi:
	rm -rf api

clean-local: cleanapi cleancoverage
	find . -name "*~" -delete
	rm -rf *.stamp .*.stamp compile scriptmacros.m4 $(BINDIR)
