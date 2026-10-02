/**
 * @file DiracAngle.h
 *
 * @brief Angular machinery shared by the retarded (multipole) Dirac
 *        transition codes
 *
 * The spinor-harmonic decomposition and the plane-wave angular integral
 * used by the retarded matrix-element assembly in DiracRate.cxx (bound-bound
 * rates) and DiracRR.cxx (bound-free radiative recombination).  The
 * conventions here are the single source of truth for the phase factors of
 * the retarded operator e^{-i k r cosT}; both translation units must agree
 * exactly.
 *
 * The implementations live in DiracAngle.cxx.
 *
 * $Id: DiracAngle.h 2122 2026-09-28 11:42:05Z iamp $
 // SPDX-License-Identifier: MIT
 */

#pragma once

#include <complex>
#include <vector>

/**
 * One component of a spinor harmonic chi_{kappa,m}: (coeff, l, m).
 */
struct AngComp {
  double coeff;
  int l;
  int m;
};

/**
 * @brief the two components of chi_{kappa,m}
 *
 * chi_{kappa,m} = CG(l, m-1/2, 1/2, 1/2; j, m) Y_{l,m-1/2} chi_up
 *               + CG(l, m+1/2, 1/2, -1/2; j, m) Y_{l,m+1/2} chi_down
 *
 * with l = l(kappa), j = j(kappa) and m = ms2/2.  A component with
 * |m| > l does not exist; its coefficient is then set to zero.
 *
 * @param kappa relativistic angular quantum number
 * @param ms2   2*m (an odd integer)
 * @param c1    (output) first (spin-up) spinor component
 * @param c2    (output) second (spin-down) spinor component
 */
void spinor_harm(int kappa, int ms2, AngComp &c1, AngComp &c2);

/**
 * @brief coefficients ac[L] such that
 *        int dOm Y*_{l1,m1} Y_{l2,m2} e^{-is cosT} = sum_L ac[L] j_L(s)
 *
 * With e^{-i k r cos theta} = sum_L (2L+1) (-i)^L j_L(k r) P_L(cos theta),
 *   ac[L] = (-i)^L (-1)^{m1} sqrt((2l1+1)(2l2+1))
 *           <l1 0, l2 0 | L 0> <l1 -m1, l2 m2 | L 0>.
 * ac[L] is zero unless m1 == m2 and |l1-l2| <= L <= l1+l2 with
 * l1+l2+L even.
 *
 * @param l1 orbital angular momentum of the first spherical harmonic
 * @param m1 magnetic quantum number of the first spherical harmonic
 * @param l2 orbital angular momentum of the second spherical harmonic
 * @param m2 magnetic quantum number of the second spherical harmonic
 *
 * @return the vector ac[L] for L = 0..l1+l2
 */
std::vector<std::complex<double> > ang_int(int l1, int m1, int l2, int m2);

/**
 * @brief spherical Bessel function j_L(x) for x >= 0
 *
 * Three regimes: a power series for x < 0.5, upward recurrence for
 * L <= x, and downward (Miller) recurrence for L > x.  Accurate to
 * machine precision over the whole range of arguments produced by the
 * radial grids of the retarded codes.
 *
 * @param L order of the spherical Bessel function
 * @param x argument
 *
 * @return j_L(x)
 */
double sbesj(int L, double x);

