/**
 * @file hypergeometric.h
 *
 * @brief definitions of hypergeometric functions
 *
 * @par CREATION
 * @author Stefan Schippers
 * @date 2026
 *
 * @par VERSION
 * @verbatim
 * $Id: hypergeometric.h 2118 2026-09-22 09:12:57Z iamp $
 // SPDX-License-Identifier: MIT
 * @endverbatim
 */

#pragma once

#include <complex>

// FLINT (arb) arbitrary-precision ball arithmetic for the exact 1F1 and 2F1
// engines. CMake links FLINT and defines HYDROCAL_HAVE_FLINT when pkg-config
// finds it (see CMakeLists.txt), and the __has_include guard additionally
// protects against a stale define with the headers absent. The
// HYDROCAL_HAVE_ACB_HYPGEOM macro is defined here rather than only inside
// hypergeometric.cxx so that every translation unit including this header
// (notably DiracRR.cxx, which gates the fully-retarded cross section and the
// analytic Gordon element on it) sees the same availability signal.
#if defined(HYDROCAL_HAVE_FLINT) && __has_include(<flint/acb.h>)
#define HYDROCAL_HAVE_ACB_HYPGEOM 1
#endif


///////////////////////////////////////////////////////////////////////////////
/** 
 * @brief  Complex gamma function
 */
std::complex<double> cgamma_cmplx(const std::complex<double> &z);


///////////////////////////////////////////////////////////////////////////////
/** 
 * @brief  Complex log-gamma function
 */
std::complex<double> clngamma_cmplx(const std::complex<double> &z);


/////////////////////////////////////////////////////////////////
/**
 * @brief Complex confluent hypergeometric function 1F1(a; b; z)
 */
void hypconfl_cmplx(double a_re, double a_im, double b, double z_re,
		    double z_im, std::complex<double> &res);

////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Complex Gauss hypergeometric function 2F1(a,b; c; z)
 */
bool hypgeo2_cmplx(const std::complex<double> &a, const std::complex<double> &b,
		   const std::complex<double> &c, const std::complex<double> &z,
		   std::complex<double> &res);
