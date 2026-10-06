// SPDX-License-Identifier: MIT
/**
 * @file hypergeometric.cxx
 *
 * @brief coding of hypergeometric functions
 *
 * @par CREATION
 * @author Stefan Schippers
 * @date 2026
 *
 * @par VERSION
 * @verbatim
 * @endverbatim
 */

#include "hydroconst.h"
#include "hypergeometric.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

// FLINT (arb) arbitrary-precision ball arithmetic for the exact 1F1 and 2F1 engines.
// Detection mirrors the Boost pattern: CMake links the library and defines
// HYDROCAL_HAVE_FLINT when pkg-config finds it, and the __has_include guard
// here additionally protects against a stale define with the headers absent.
// HYDROCAL_HAVE_ACB_HYPGEOM itself is defined in hypergeometric.h (included
// above) so that every translation unit shares one definition of the signal.
#if defined(HYDROCAL_HAVE_FLINT) && __has_include(<flint/acb.h>)
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#endif

#if defined(__SIZEOF_FLOAT128__) && !defined(HYDROCAL_HAVE_ACB_HYPGEOM)
#define HYDROCAL_HAVE_FLOAT128 1
struct cq { // complex arithmetic for __float128
  __float128 re, im;
};
static cq cq_add(cq a, __float128 r) { return {a.re + r, a.im}; }
static cq cq_div(cq a, __float128 r) { return {a.re / r, a.im / r}; }
static cq cq_mul(cq a, cq b) {
  return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
static __float128 cq_abs2(cq a) { return a.re * a.re + a.im * a.im; }
static cq cq_divc(cq a, cq b) {
  __float128 d = b.re * b.re + b.im * b.im;
  return {(a.re * b.re + a.im * b.im) / d,
          (a.im * b.re - a.re * b.im) / d};
}
#endif

using namespace std;



///////////////////////////////////////////////////////////////////////////////
/** 
 * @brief  Complex gamma function via the Lanczos approximation
 *
 * (g=7, n=9, Godfrey coefficients).
 *
 * All arguments used here (a, b-a, a+1, b-(a+1), b, 2s, ...)
 * have positive real part, but the reflection formula is included for robustness.
 */
complex<double> cgamma_cmplx(const complex<double> &z) {
  static const double coef[9] = {
      0.99999999999980993, 676.5203681218851, -1259.1392167224028,
      771.32342877765313,  -176.61502916214059, 12.507343278686905,
      -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
  if (z.real() <= 0.0) {
    complex<double> pp = hydroconst::pi / sin(hydroconst::pi * z) /
                         cgamma_cmplx(1.0 - z);
    return pp;
  }
  complex<double> zz = z - 1.0;
  complex<double> x = coef[0];
  for (int i = 1; i < 9; i++) x += coef[i] / (zz + (double)i);
  complex<double> t = zz + 7.5;
  return sqrt(2.0 * hydroconst::pi) * pow(t, zz + 0.5) * exp(-t) * x;
}


///////////////////////////////////////////////////////////////////////////////
/** 
 * @brief  Complex log-gamma via the same Lanczos series.
 *
 * Returning the logarithm
 *
 * Keeps large-|Im z| arguments representable: |Gamma(s+iy)| falls off like
 * e^{-pi|y|/2}, which underflows double for |y| of a few hundred, whereas
 * log|Gamma| stays an ordinary O(|y|) number.  The closed-form normalisation
 * needs this because Im(a) = s b1/k ~ -pi*eta can reach thousands at low
 * energy, and the branch amplitudes C_in, C_out stay O(1) only if the
 * e^{-pi|y|/2} falloff of |Gamma(a)| is combined with the compensating
 * e^{+pi|y|/2} of |(-2ik)^{a-2s}| *inside the logarithm* rather than formed
 * separately (where both factors would flush to zero).
 */
complex<double> clngamma_cmplx(const complex<double> &z) {
  static const double coef[9] = {
      0.99999999999980993, 676.5203681218851, -1259.1392167224028,
      771.32342877765313,  -176.61502916214059, 12.507343278686905,
      -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
  if (z.real() <= 0.0) {
    return log(hydroconst::pi) - log(sin(hydroconst::pi * z)) -
           clngamma_cmplx(1.0 - z);
  }
  complex<double> zz = z - 1.0;
  complex<double> x = coef[0];
  for (int i = 1; i < 9; i++) x += coef[i] / (zz + (double)i);
  complex<double> t = zz + 7.5;
  return 0.5 * log(2.0 * hydroconst::pi) + (zz + 0.5) * log(t) - t + log(x);
}

///////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Confluent hypergeometric function 1F1 by direct Taylor series in double precision
 *
 * @param a_re real part of complex a
 * @param a_im imaginary part of complex a
 * @param b real
 * @param z_re real part of complex z
 * @param z_rim imaginary part of complex z
 * @param res on exit: complex value of 1F1(a; b; z) 
 *
 * @return true if return value is trusted
 * @return false if return value is not trusted
 *
 * 1F1(a; b; z) = sum_n (a)_n z^n / ((b)_n n!)
 *
 * The series is trusted only when its roundoff estimate
 * peak*eps_double/|1F1| (peak = largest |term| encountered, a proxy for the
 * accumulated rounding error) is <= 1e-8; for large |z| the intermediate
 * terms grow ~e^|z| before decaying and the sum is lost to cancellation, in
 * which case this returns false and the caller falls through to the next
 * stage.  Returns false immediately (fast) if the terms overflow.
 */
static bool hypconfl_series_double(double a_re, double a_im, double b,
                                   double z_re, double z_im,
                                   complex<double> &res) {
  const int NMAX = 200000;
  const double eps = numeric_limits<double>::epsilon();
  complex<double> a(a_re, a_im), z(z_re, z_im);
  complex<double> sum(1.0, 0.0), term(1.0, 0.0);
  double peak = 1.0;
  bool conv = false;
  for (int n = 0; n < NMAX; n++) {
    term *= ((a + (double)n) / (b + (double)n)) * (z / (double)(n + 1));
    sum += term;
    double m2 = norm(term);
    if (m2 > peak * peak) peak = sqrt(m2);
    double s2 = norm(sum);
    if (!(s2 > 0.0) || !(m2 <= 1e308)) return false; // nan/inf -> no trust
    if (n > 0 && m2 <= 64.0 * eps * eps * s2) { conv = true; break; }
  }
  if (!conv) return false;
  if (peak * eps / abs(sum) > 1e-8) return false;
  res = sum;
  return true;
}

//////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Confluent hypergeometric function 1F1 by direct Taylor series in quad precision
 *
 * @param a_re real part of complex a
 * @param a_im imaginary part of complex a
 * @param b real
 * @param z_re real part of complex z
 * @param z_rim imaginary part of complex z
 * @param res on exit: complex value of 1F1(a; b; z) 
 *
 * @return true if return value is trusted
 * @return false if return value is not trusted
 *
 * The same Taylor series in quad precision (mantissa rounding 2^-110) and
 * the Poincare asymptotic expansion. Trusted where the roundoff estimate
 * peak*eps_quad/|1F1| <= 1e-10; beyond that (very large |z|, terms overflowing even quad)
 * returns false.
 *
 * Superseeded by FLINT when available.
 */
#ifdef HYDROCAL_HAVE_FLOAT128
static bool hypconfl_series_quad(double a_re, double a_im, double b,
                                 double z_re, double z_im,
                                 complex<double> &res) {
  const int NMAX = 200000;
  const __float128 epsq = (__float128)7.6e-34L; // ~2^-110
  cq a{(__float128)a_re, (__float128)a_im};
  const __float128 bb = (__float128)b;
  const cq z{(__float128)z_re, (__float128)z_im};
  cq sum{1, 0}, term{1, 0};
  __float128 peak2 = 1.0;
  bool conv = false;
  for (int n = 0; n < NMAX; n++) {
    __float128 nf = (__float128)n;
    cq coef = cq_div(cq_add(a, nf), bb + nf);
    cq dz = cq_div(z, nf + 1.0);
    term = cq_mul(cq_mul(term, coef), dz);
    sum.re += term.re;
    sum.im += term.im;
    __float128 m2 = cq_abs2(term);
    if (m2 > peak2) peak2 = m2;
    __float128 s2 = cq_abs2(sum);
    if (!(s2 <= (__float128)1e4900L)) return false; // overflow -> no trust
    if (s2 != s2) return false;                     // nan
    if (n > 0 && m2 <= (__float128)1e-60L * s2) { conv = true; break; }
  }
  if (!conv) return false;
  // roundoff estimate peak*epsq/|sum| <= 1e-10, in squared form
  if (peak2 * epsq * epsq > (__float128)1e-20L * cq_abs2(sum)) return false;
  res = complex<double>((double)sum.re, (double)sum.im);
  return true;
}
#endif

//////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Confluent hypergeometric function 1F1 by Poincare asymptotic expansion (DLMF 13.7.1)
 *
 * @param a_re real part of complex a
 * @param a_im imaginary part of complex a
 * @param b real
 * @param z_re real part of complex z
 * @param z_rim imaginary part of complex z
 * @param res on exit: complex value of 1F1(a; b; z) 
 *
 * truncated at its smallest term. Used in the far field (large |z|) where it is valid and
 * the Taylor series has no chance.
 *
 * Superseeded by FLINT when available.
 */
#if !defined(HYDROCAL_HAVE_ACB_HYPGEOM)
static void hypconfl_asym(double a_re, double a_im, double b, double z_re,
                          double z_im, complex<double> &res) {
  const int KMAX = 100000;
  complex<double> a(a_re, a_im), z(z_re, z_im);
  complex<double> ba = b - a;
  complex<double> gb = clngamma_cmplx(b);
  complex<double> gba = clngamma_cmplx(ba);
  complex<double> ga = clngamma_cmplx(a);

  complex<double> s1(1.0, 0.0), t1(1.0, 0.0);
  double mprev = 1.0;
  for (int k = 1; k < KMAX; k++) {
    t1 *= ((a + (double)(k - 1)) * (a - b + (double)k)) /
          ((double)k * (-z));
    double m = abs(t1);
    if (m > mprev) break;
    s1 += t1;
    mprev = m;
  }

  complex<double> s2(1.0, 0.0), t2(1.0, 0.0);
  mprev = 1.0;
  for (int k = 1; k < KMAX; k++) {
    t2 *= ((ba + (double)(k - 1)) * ((double)k - a)) / ((double)k * z);
    double m = abs(t2);
    if (m > mprev) break;
    s2 += t2;
    mprev = m;
  }

  res = exp(gb - gba) * pow(-z, -a) * s1 +
        exp(gb - ga) * exp(z) * pow(z, a - b) * s2;
}
#endif

/////////////////////////////////////////////////////////////////
/**
 * @brief Confluent hypergeometric function 1F1 from FLINT
 *
 * @param a_re real part of complex a
 * @param a_im imaginary part of complex a
 * @param b real
 * @param z_re real part of complex z
 * @param z_rim imaginary part of complex z
 * @param res on exit: complex value of 1F1(a; b; z) 
 *
 * exact ball arithmetic for the deep-cancellation regime.
 */
#ifdef HYDROCAL_HAVE_ACB_HYPGEOM
static void hypconfl_flint(double a_re, double a_im, double b, double z_re,
                           double z_im, complex<double> &res) {
  acb_t aa, bb, zz, rres;
  acb_init(aa);
  acb_init(bb);
  acb_init(zz);
  acb_init(rres);
  acb_set_d_d(aa, a_re, a_im);
  acb_set_d(bb, b);
  acb_set_d_d(zz, z_re, z_im);
  // 100 bits (~30 digits) is far beyond the ~1e-9 relative accuracy the
  // closed-form wavefunction needs, without the cost of a much higher
  // precision ball computation on a large radial grid.
  acb_hypgeom_1f1(rres, aa, bb, zz, 0, 100);
  res.real(arf_get_d(arb_midref(acb_realref(rres)), ARF_RND_NEAR));
  res.imag(arf_get_d(arb_midref(acb_imagref(rres)), ARF_RND_NEAR));
  acb_clear(aa);
  acb_clear(bb);
  acb_clear(zz);
  acb_clear(rres);
}
#endif

/////////////////////////////////////////////////////////////////
/**
 * @brief Confluent hypergeometric function 1F1 (wrapper)
 *
 * @param a_re real part of complex a
 * @param a_im imaginary part of complex a
 * @param b real
 * @param z_re real part of complex z
 * @param z_rim imaginary part of complex z
 * @param res on exit: complex value of 1F1(a; b; z) 
 *
 * Engine-selecting wrapper.  he cheap direct double Taylor series is always
 * tried first (it is accurate for the vast majority of grid points and costs
 * microseconds); the deep-cancellation points then go to FLINT when
 * available, else to the quad/asymptotic stages.
 */
void hypconfl_cmplx(double a_re, double a_im, double b, double z_re,
                           double z_im, complex<double> &res) {
  if (hypconfl_series_double(a_re, a_im, b, z_re, z_im, res)) return;
#ifdef HYDROCAL_HAVE_ACB_HYPGEOM
  hypconfl_flint(a_re, a_im, b, z_re, z_im, res);
#else
#ifdef HYDROCAL_HAVE_FLOAT128
  if (hypconfl_series_quad(a_re, a_im, b, z_re, z_im, res)) return;
#endif
  hypconfl_asym(a_re, a_im, b, z_re, z_im, res);
#endif
}

////////////////////////////////////////////////////////////////////////////
/**
 * @brief Gauss hypergeometric function 2F1(a, b; c; z) for complex a, b, c, z.
 *
 * @param a complex
 * @param b complex
 * @param c complex
 * @param z complex
 * @param res on exit compex value of 2F1(a,b; c; z)
 *
 * @return true 
 *
 * The analytic (relativistic-Gordon) bound-free radial integral
 * (dirac_rr_bf_analytic below) needs 2F1(a, nu; 2s; z) and
 * 2F1(a+1, nu; 2s+1; z) with a = s(1+ib1/k), real nu, and the argument
 *   z = -2ik/(lambda - ik),  |z| = 2k/sqrt(lambda^2 + k^2) < 2,
 * which leaves the unit disk of the Gauss series at high energy
 * (k >> lambda) -- precisely the deep-cancellation regime the analytic
 * element is meant to fix.  Two engines:
 *   - FLINT (arb ball arithmetic, HYDROCAL_HAVE_ACB_HYPGEOM): exact, with
 *     automatic analytic continuation; the default when available.
 *   - a portable double/quad fallback (no FLINT): Gauss series with the
 *     same roundoff estimate as the 1F1 engine, the Pfaff transformation
 *     to z/(z-1), and the reciprocal connection formula (DLMF 15.8.2) for
 *     the large-|z| branch.  Returns false when no stage is trustworthy,
 *     so the caller can fall back to the numerical path.
 */

////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Gauss hypergeometric function 2F1 for complex arguments by Gauss series in double precision 
 *
 * @param a complex
 * @param b complex
 * @param c complex
 * @param z complex
 * @param res on exit compex value of 2F1(a,b; c; z)
 *
 * @return true if return value is trusted
 * @return false if return value is not trusted
 *
 * trusted only when the roundoff estimate peak*eps/|2F1| <= 1e-8, mirroring hypconfl_series_double.
 */
static bool hypgeo2_gauss_double(const complex<double> &a,
				 const complex<double> &b,
				 const complex<double> &c,
				 const complex<double> &z,
				 complex<double> &res) {
  const int NMAX = 200000;
  const double eps = numeric_limits<double>::epsilon();
  complex<double> sum(1.0, 0.0), term(1.0, 0.0);
  double peak = 1.0;
  bool conv = false;
  for (int n = 0; n < NMAX; n++) {
    term *= ((a + (double)n) * (b + (double)n) /
             ((c + (double)n) * (double)(n + 1))) * z;
    sum += term;
    double m2 = norm(term);
    if (m2 > peak * peak) peak = sqrt(m2);
    double s2 = norm(sum);
    if (!(s2 > 0.0) || !(m2 <= 1e308)) return false;  // nan/inf -> no trust
    if (n > 0 && m2 <= 64.0 * eps * eps * s2) { conv = true; break; }
  }
  if (!conv) return false;
  if (peak * eps / abs(sum) > 1e-8) return false;
  res = sum;
  return true;
}

#ifdef HYDROCAL_HAVE_FLOAT128
////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Gauss hypergeometric function 2F1 for complex arguments by Gauss series in quad precision 
 *
 * @param a complex
 * @param b complex
 * @param c complex
 * @param z complex
 * @param res on exit compex value of 2F1(a,b; c; z)
 *
 * @return true if return value is trusted
 * @return false if return value is not trusted
 *
 * Gauss series in quad precision, for the band around |z| = 1 where the
 * double series passes its convergence test only with an unacceptable
 * roundoff estimate.  Trusted where peak*eps_quad/|2F1| <= 1e-10.
 */
static bool hypgeo2_gauss_quad(const complex<double> &a,
                               const complex<double> &b,
                               const complex<double> &c,
                               const complex<double> &z,
                               complex<double> &res) {
  const int NMAX = 200000;
  const __float128 epsq = (__float128)7.6e-34L;
  const cq ca{(__float128)a.real(), (__float128)a.imag()};
  const cq cb{(__float128)b.real(), (__float128)b.imag()};
  const cq cc{(__float128)c.real(), (__float128)c.imag()};
  const cq cz{(__float128)z.real(), (__float128)z.imag()};
  cq sum{1, 0}, term{1, 0};
  __float128 peak2 = 1.0;
  bool conv = false;
  for (int n = 0; n < NMAX; n++) {
    __float128 nf = (__float128)n;
    cq num = cq_mul(cq_add(ca, nf), cq_add(cb, nf));
    cq den = cq_mul(cq_add(cc, nf), cq{nf + 1.0, 0});
    term = cq_mul(cq_divc(num, den), cz);
    sum.re += term.re;
    sum.im += term.im;
    __float128 m2 = cq_abs2(term);
    if (m2 > peak2) peak2 = m2;
    __float128 s2 = cq_abs2(sum);
    if (!(s2 <= (__float128)1e4900L)) return false;  // overflow -> no trust
    if (s2 != s2) return false;                      // nan
    if (n > 0 && m2 <= (__float128)1e-60L * s2) { conv = true; break; }
  }
  if (!conv) return false;
  if (peak2 * epsq * epsq > (__float128)1e-20L * cq_abs2(sum)) return false;
  res = complex<double>((double)sum.re, (double)sum.im);
  return true;
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Gauss hypergeometric function 2F1 for complex arguments inside unit disk
 *
 * @param a complex
 * @param b complex
 * @param c complex
 * @param z complex
 * @param res on exit compex value of 2F1(a,b; c; z)
 *
 * @return true if return value is trusted
 * @return false if return value is not trusted
 *
 * Portable 2F1: Gauss series inside the unit disk, Pfaff transformations
 * (DLMF 15.8.4) and the reciprocal connection formula (DLMF 15.8.2, valid
 * for |arg(-z)| < pi, true for every z this code passes here). Depth is
 * bounded so the recursive transform calls always terminate.
 */
static bool hypgeo2_portable(const complex<double> &a,
                             const complex<double> &b,
                             const complex<double> &c,
                             const complex<double> &z,
                             complex<double> &res, int depth) {
  if (depth > 3) return false;
  if (abs(z) <= 0.95) {
    if (hypgeo2_gauss_double(a, b, c, z, res)) return true;
#ifdef HYDROCAL_HAVE_FLOAT128
    if (hypgeo2_gauss_quad(a, b, c, z, res)) return true;
#endif
    return false;
  }
  {
    complex<double> w = z / (z - 1.0);
    if (abs(w) < 1.0) {
      complex<double> t;
      if (hypgeo2_portable(a, c - b, c, w, t, depth + 1)) {
        res = pow(1.0 - z, -a) * t;
        return true;
      }
      if (hypgeo2_portable(c - a, b, c, w, t, depth + 1)) {
        res = pow(1.0 - z, -b) * t;
        return true;
      }
    }
  }
  if (abs(z) >= 1.2) {
    complex<double> u = 1.0 / z;
    complex<double> uz = -z;
    complex<double> A1 = cgamma_cmplx(c) * cgamma_cmplx(b - a) /
                         (cgamma_cmplx(b) * cgamma_cmplx(c - a));
    complex<double> A2 = cgamma_cmplx(c) * cgamma_cmplx(a - b) /
                         (cgamma_cmplx(a) * cgamma_cmplx(c - b));
    complex<double> t1, t2;
    if (hypgeo2_portable(a, a + 1.0 - c, a + 1.0 - b, u, t1, depth + 1) &&
        hypgeo2_portable(b, b + 1.0 - c, b + 1.0 - a, u, t2, depth + 1)) {
      res = A1 * exp(-a * log(uz)) * t1 + A2 * exp(-b * log(uz)) * t2;
      return true;
    }
  }
  return false;
}

#ifdef HYDROCAL_HAVE_ACB_HYPGEOM
////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Gauss hypergeometric function 2F1 from FLINT
 *
 * @param a complex
 * @param b complex
 * @param c complex
 * @param z complex
 * @param res on exit compex value of 2F1(a,b; c; z)
 *
 * FLINT engine: exact ball arithmetic with automatic continuation.
 */
static void hypgeo2_flint(const complex<double> &a, const complex<double> &b,
                          const complex<double> &c, const complex<double> &z,
                          complex<double> &res) {
  acb_t aa, bb, cc, zz, rres;
  acb_init(aa);
  acb_init(bb);
  acb_init(cc);
  acb_init(zz);
  acb_init(rres);
  acb_set_d_d(aa, a.real(), a.imag());
  acb_set_d_d(bb, b.real(), b.imag());
  acb_set_d_d(cc, c.real(), c.imag());
  acb_set_d_d(zz, z.real(), z.imag());
  acb_hypgeom_2f1(rres, aa, bb, cc, zz, 0, 256);
  res.real(arf_get_d(arb_midref(acb_realref(rres)), ARF_RND_NEAR));
  res.imag(arf_get_d(arb_midref(acb_imagref(rres)), ARF_RND_NEAR));
  acb_clear(aa);
  acb_clear(bb);
  acb_clear(cc);
  acb_clear(zz);
  acb_clear(rres);
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Gauss hypergeometric function 2F1 (wrapper)
 *
 * @param a complex
 * @param b complex
 * @param c complex
 * @param z complex
 * @param res on exit compex value of 2F1(a,b; c; z)
 *
 * Engine-selecting wrapper. Returns false only if the portable engine was
 * unable to produce a trustworthy value (FLINT always succeeds).
 */
bool hypgeo2_cmplx(const complex<double> &a, const complex<double> &b,
		   const complex<double> &c, const complex<double> &z,
		   complex<double> &res) {
#ifdef HYDROCAL_HAVE_ACB_HYPGEOM
  hypgeo2_flint(a, b, c, z, res);
  return true;
#else
  return hypgeo2_portable(a, b, c, z, res, 0);
#endif
}

