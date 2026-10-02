/**
 * @file DiracAngle.cxx
 *
 * @brief Implementation of the angular machinery shared by the fully
 *        retarded bound-bound and bound-free transition codes
 *
 * @par CREATION
 * @author Stefan Schippers
 * @date 2026
 *
 * @par VERSION
 * @verbatim
 * $Id: DiracAngle.cxx 2112 2026-08-14 17:40:06Z iamp $
 // SPDX-License-Identifier: MIT
 * @endverbatim
 */
#include "DiracAngle.h"
#include "clebsch.h"
#include <cmath>

using namespace std;

//////////////////////////////////////////////////////////////////////////
/**
 * @brief decompose the spinor spherical harmonic chi_{kappa,m} into its
 *        spin-up and spin-down components
 *
 * chi_{kappa,m} = CG(l, m-1/2, 1/2, 1/2; j, m) Y_{l,m-1/2} chi_up
 *               + CG(l, m+1/2, 1/2, -1/2; j, m) Y_{l,m+1/2} chi_down
 *
 * with l = l(kappa), j = j(kappa) and m = ms2/2.  The returned m quantum
 * numbers (ms2 +- 1)/2 are always integers.
 */
void spinor_harm(int kappa, int ms2, AngComp &upper, AngComp &lower) {
  double j = fabs((double)kappa) - 0.5;
  int l = (kappa > 0) ? kappa : -kappa - 1;

  if (abs(ms2) > (int)(2.0 * j)) {
    upper.coeff = 0.0;
    upper.l = l;
    upper.m = 0;
    lower.coeff = 0.0;
    lower.l = l;
    lower.m = 0;
    return;
  }

  int upper_m = (ms2 - 1) / 2; // m - 1/2 (integer)
  int lower_m = (ms2 + 1) / 2; // m + 1/2 (integer)
  double m = 0.5 * ms2;

  upper.coeff =
      (abs(upper_m) <= l) ? CG(l, 0.5, upper_m, 0.5, j, m) : 0.0;
  upper.l = l;
  upper.m = upper_m;

  lower.coeff =
      (abs(lower_m) <= l) ? CG(l, 0.5, lower_m, -0.5, j, m) : 0.0;
  lower.l = l;
  lower.m = lower_m;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief angular integral of two spherical harmonics against the retarded
 *        plane wave, expanded in spherical Bessel functions
 *
 * Using e^{-i k r cos theta} = sum_L (2L+1) (-i)^L j_L(k r) P_L(cos theta)
 * and P_L(cos theta) = sqrt(4 pi/(2L+1)) Y_{L0}(theta, phi),
 *   int dOmega Y*_{l1 m1} Y_{l2 m2} e^{-i k r cos theta}
 *     = sum_L c_L j_L(k r),
 *   c_L = (-i)^L (-1)^{m1} sqrt((2l1+1)(2l2+1))
 *         <l1 0, l2 0 | L 0> <l1 -m1, l2 m2 | L 0>,
 * where the (2L+1) of the plane-wave expansion and the 1/sqrt(2L+1) of
 * P_L = sqrt(4 pi/(2L+1)) Y_{L0} cancel against the two sqrt(2L+1) of the
 * Clebsch-Gordan coefficients.
 * The second Clebsch-Gordan coefficient vanishes unless m1 == m2.
 */
vector<complex<double>> ang_int(int l1, int m1, int l2, int m2) {
  vector<complex<double>> c(l1 + l2 + 1, complex<double>(0.0, 0.0));

  // Out-of-range magnetic quantum numbers (possible for the vanishing
  // components of a spinor harmonic, e.g. l = 0 with m = -1) make the
  // Clebsch-Gordan machinery produce an invalid 3j symbol; the integral
  // is zero in that case.
  if (abs(m1) > l1 || abs(m2) > l2) return c;

  double pref = sqrt((2.0 * l1 + 1.0) * (2.0 * l2 + 1.0));
  int Lmin = abs(l1 - l2);
  for (int L = Lmin; L <= l1 + l2; L++) {
    if ((l1 + l2 + L) % 2 != 0) continue; // <l1 0, l2 0 | L 0> = 0
    double g00 = CG(l1, l2, 0, 0, L, 0);
    if (g00 == 0.0) continue;
    double gmm = CG(l1, l2, -m1, m2, L, 0);
    if (gmm == 0.0) continue;

    double amp = pref * g00 * gmm;
    if (m1 & 1) amp = -amp; // (-1)^m1

    // (-i)^L:  1, -i, -1, i, 1, ...
    complex<double> ipow(1.0, 0.0);
    if (L % 4 == 1)
      ipow = complex<double>(0.0, -1.0);
    else if (L % 4 == 2)
      ipow = complex<double>(-1.0, 0.0);
    else if (L % 4 == 3)
      ipow = complex<double>(0.0, 1.0);

    c[L] = amp * ipow;
  }
  return c;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief spherical Bessel function j_L(x) for x >= 0
 *
 * Three regimes:
 *   - x < 0.5:    power series (rapid, well-conditioned)
 *   - L <= x:     upward recurrence from exact j_0, j_1 (stable for L < x)
 *   - L >  x:     downward (Miller) recurrence from N = L + 30
 */
double sbesj(int L, double x) {
  const double pi = 3.14159265358979323846;
  if (L < 0) return 0.0;
  if (x < 0.0) x = -x;

  if (x < 0.5) {
    // j_L(x) = (sqrt(pi)/2) (x/2)^L sum_k (-1)^k (x/2)^(2k)
    //                                      / (k! Gamma(L+k+3/2))
    double hx2 = 0.25 * x * x;
    double j = 0.5 * sqrt(pi) * pow(0.5 * x, (double)L) /
               tgamma(L + 1.5);
    double term = 1.0, sum = 1.0;
    for (int k = 1; k <= 2000; k++) {
      term *= -hx2 / (k * (L + k + 0.5));
      sum += term;
      if (fabs(term) < fabs(sum) * 1e-16) break;
    }
    return j * sum;
  }

  if (L <= x) {
    // Upward recurrence j_{l+1} = (2l+1)/x j_l - j_{l-1}
    double j0 = sin(x) / x;
    if (L == 0) return j0;
    double j1 = (sin(x) - x * cos(x)) / (x * x);
    if (L == 1) return j1;
    double jm1 = j0, jm = j1;
    for (int l = 1; l < L; l++) {
      double jp = (2.0 * l + 1.0) / x * jm - jm1;
      jm1 = jm;
      jm = jp;
    }
    return jm;
  }

  // Downward (Miller) recurrence, stable in the decaying regime L > x.
  // Seed j[N+1] = 0, j[N] = 1e-300; the common scale is fixed by
  // normalising the computed j_0 against sin(x)/x.
  const int N = L + 30;
  double jNp1 = 0.0, jN = 1e-300;
  double jL = 0.0;
  for (int n = N; n >= 1; n--) {
    double jnm1 = (2.0 * n + 1.0) / x * jN - jNp1;
    jNp1 = jN;
    jN = jnm1;
    if (n - 1 == L) jL = jN;
  }
  return jL * (sin(x) / x / jN);
}
