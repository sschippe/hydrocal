// SPDX-License-Identifier: MIT
/**
 * @file DiracRate.cxx
 *
 * @brief Relativistic (Dirac) hydrogenic E1, E2, M1 transition rates,
 *        lifetimes and branching ratios
 *
 * Uses analytic Dirac hydrogenic wavefunctions expressed in terms
 * of confluent hypergeometric functions.  The radial multipole
 * integrals and the wavefunction normalisation are computed
 * analytically via polynomial expansion & Gamma-function
 * integrals.
 *
 * @par CREATION
 * @author Stefan Schippers
 * @date 2026
 *
 * @par VERSION
 * @verbatim
 * @endverbatim
 */
#include "DiracRate.h"
#include "DiracAngle.h"
#include "clebsch.h"
#include "hydroconst.h"
#include "hydromath.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>


using namespace std;

//////////////////////////////////////////////////////////////////////////
/**
 * @brief converts orbital and total angular momentum to the Dirac quantum
 *        number kappa
 *
 * @param l orbital angular momentum quantum number
 * @param j total angular momentum (half-integer, e.g. 1.5)
 *
 * @return kappa = -(l+1) for j = l+1/2, kappa = +l for j = l-1/2
 */
int kappa_from_lj(int l, double j) {
  if (l == 0) return -1;
  if (j > l - 0.25)
    return -(l + 1);
  else
    return l;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief converts the Dirac quantum number kappa to l and j
 *
 * @param kappa relativistic angular quantum number
 * @param l     (output) orbital angular momentum quantum number
 * @param j     (output) total angular momentum (half-integer, e.g. 1.5)
 */
void lj_from_kappa(int kappa, int &l, double &j) {
  if (kappa < 0) {
    l = -kappa - 1;
    j = l + 0.5;
  } else {
    l = kappa;
    j = l - 0.5;
  }
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief checks whether a state (n, kappa) exists
 *
 * @param n     principal quantum number
 * @param kappa relativistic angular quantum number
 *
 * @return 1 if l < n, i.e. the state is physical, 0 otherwise
 */
bool valid_kappa(int n, int kappa) {
  int l;
  double j;
  lj_from_kappa(kappa, l, j);
  return l < n;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Dirac (relativistic) binding energy
 *
 * @param z     nuclear charge
 * @param n     principal quantum number
 * @param kappa relativistic angular quantum number
 *
 * @return binding energy in eV (negative)
 */
double dirac_binding_energy(double z, int n, int kappa) {
  using hydroconst::alpha;
  using hydroconst::mec2_eV;
  double za = z * alpha;
  int nr = n - abs(kappa);
  double gamma = sqrt(kappa * kappa - za * za);
  double t = za / (nr + gamma);
  return mec2_eV * (1.0 / sqrt(1.0 + t * t) - 1.0);
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief number of Dirac states up to principal quantum number n
 *
 * @param n principal quantum number
 *
 * @return n*n
 */
int ndirac_states(int n) { return n * n; }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief index of the state (n, kappa) in the counting scheme of all
 *        Dirac states up to n
 *
 * @param n     principal quantum number
 * @param kappa relativistic angular quantum number
 *
 * @return sequential index (0-based), inverse of kappa_from_index
 */
int dirac_state_index(int n, int kappa) {
  int prev = (n - 1) * (n - 1);
  int l;
  double j;
  lj_from_kappa(kappa, l, j);
  int count = 0;
  for (int lp = 0; lp < l; lp++) count += (lp == 0) ? 1 : 2;
  if (l > 0 && kappa > 0) count++;
  return prev + count;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief inverse mapping of dirac_state_index: state index to (n, kappa)
 *
 * @param idx    sequential state index (0-based)
 * @param n      (output) principal quantum number
 * @param kappa  (output) relativistic angular quantum number
 */
void kappa_from_index(int idx, int &n, int &kappa) {
  n = (int)floor(sqrt((double)idx)) + 1;
  int offset = (n - 1) * (n - 1);
  int k = idx - offset;
  for (int l = 0; l < n; l++) {
    if (l == 0) {
      if (k == 0) { kappa = -1; return; }
      k--;
    } else {
      if (k == 0) { kappa = l; return; }
      if (k == 1) { kappa = -(l + 1); return; }
      k -= 2;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief builds the polynomial coefficients of F(r) and G(r) for a Dirac
 *        state
 *
 * F(r) = sum_{k=0}^{nr}  a[k] * r^k
 * G(r) = sum_{k=0}^{nr}  b[k] * r^k
 * where nr = n - |kappa|.
 *
 * The wavefunctions are
 *   P(r) = C * sqrt(1+eps) * (2*lam*r)^gamma * exp(-lam*r) * F(r)
 *   Q(r) = C * sqrt(1-eps) * (2*lam*r)^gamma * exp(-lam*r) * G(r)
 * The polynomials are obtained from the coupled Coulomb-Dirac radial
 * equations.  In terms of u = 2*lam*r they read (see e.g. the standard
 * treatment of the hydrogenic Dirac atom)
 *   u F' + (gamma+kappa-u/2) F = (u/2 + A) G
 *   u G' + (gamma-kappa-u/2) G = (u/2 - B) F
 * with A = Z*alpha*sqrt((1-eps)/(1+eps)),
 *     B = Z*alpha*sqrt((1+eps)/(1-eps)),   A*B = (Z*alpha)^2.
 * Inserting power series F = sum a_k u^k, G = sum b_k u^k yields the
 * two-term recurrence
 *   a_k = (a_{k-1}+b_{k-1}) (k+gamma-kappa+A) / (2 k (k+2 gamma))
 *   b_k = (a_{k-1}+b_{k-1}) (k+gamma+kappa-B) / (2 k (k+2 gamma))
 * with a_0 = A, b_0 = gamma+kappa.  The series terminates at k = nr
 * because a_k + b_k ~ (k - nr), using B - A = 2 (nr+gamma) = 2 N'.
 * The returned coefficients are expanded in powers of r, i.e. they
 * carry a factor (2*lam)^k.
 *
 * @param kappa      relativistic angular quantum number
 * @param nr         n - |kappa|
 * @param gamma      sqrt(kappa^2 - (Z*alpha)^2)
 * @param eps        positive-energy amplitude factor
 * @param two_lambda 2*lam, the natural radial scale of the bound state
 * @param z          nuclear charge
 * @param a          (output) coefficients of F(r)
 * @param b_coeff    (output) coefficients of G(r)
 */
void build_FG_coeffs(int kappa, int nr, double gamma, double eps,
                            double two_lambda, double z,
                            vector<double> &a, vector<double> &b_coeff) {
  using hydroconst::alpha;
  a.assign(nr + 1, 0.0);
  b_coeff.assign(nr + 1, 0.0);

  double za = z * alpha;
  double A = za * sqrt((1.0 - eps) / (1.0 + eps));
  double B = za * sqrt((1.0 + eps) / (1.0 - eps));

  vector<double> au(nr + 1, 0.0), bu(nr + 1, 0.0);
  au[0] = A;
  bu[0] = gamma + kappa;
  for (int k = 1; k <= nr; k++) {
    double S = au[k - 1] + bu[k - 1];
    double det = 2.0 * k * (k + 2.0 * gamma);
    au[k] = S * (k + gamma - kappa + A) / det;
    bu[k] = S * (k + gamma + kappa - B) / det;
  }

  double tl = 1.0;
  for (int k = 0; k <= nr; k++) {
    if (k > 0) tl *= two_lambda;
    a[k] = au[k] * tl;
    b_coeff[k] = bu[k] * tl;
  }
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief normalisation constant of the Dirac radial functions
 *
 *   P(r) = C * sqrt(1+eps) * (2*lam*r)^gamma * exp(-lam*r) * F(r)
 *   Q(r) = C * sqrt(1-eps) * (2*lam*r)^gamma * exp(-lam*r) * G(r)
 *
 * The normalisation integral int_0^inf [P^2 + Q^2] dr = 1 fixes C.
 * It is computed analytically using the polynomial expansion.
 *
 * @param z     nuclear charge
 * @param n     principal quantum number
 * @param kappa relativistic angular quantum number
 *
 * @return the normalisation constant C
 */
double dirac_normalisation(double z, int n, int kappa) {
  using hydroconst::alpha;

  double za = z * alpha;
  int nr = n - abs(kappa);
  double gamma = sqrt(kappa * kappa - za * za);
  double Nprime = nr + gamma;
  double eps = 1.0 / sqrt(1.0 + za * za / (Nprime * Nprime));
  double lam = z / sqrt(Nprime * Nprime + za * za);   // a.u. (= sqrt(1-eps^2)/alpha)
  double two_lam = 2.0 * lam;

  vector<double> a_coef, b_coef;
  build_FG_coeffs(kappa, nr, gamma, eps, two_lam, z, a_coef, b_coef);

  // Compute  S = int_0^inf (2*lam*r)^(2*gamma) * exp(-2*lam*r)
  //                * [(1+eps)*F(r)^2 + (1-eps)*G(r)^2] dr
  //            = sum_{i,j} [(1+eps)*a_i*a_j + (1-eps)*b_i*b_j]
  //              * Gamma(2*gamma+i+j+1) / (2*lam)^(i+j+1)
  double S = 0.0;
  for (int i = 0; i <= nr; i++) {
    for (int j = 0; j <= nr; j++) {
      double comb = (1.0 + eps) * a_coef[i] * a_coef[j] +
                    (1.0 - eps) * b_coef[i] * b_coef[j];
      if (fabs(comb) < 1e-100) continue;
      double p = 2.0 * gamma + i + j + 1.0;
      S += comb * tgamma(p) / pow(two_lam, i + j + 1.0);
    }
  }

  // C = 1 / sqrt(S)  (since int[P^2+Q^2] = C^2 * S = 1)
  return 1.0 / sqrt(S);
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief fully retarded plane-wave transition rate in s^-1
 *
 * The long-wavelength multipole formulas (e.g. the dipole-only E1 formula
 *   A = (4/3) alpha^3 omega^3 (2j2+1) <kappa2||C1||kappa1>^2 R^2,
 *   R = int (P1 P2 + Q1 Q2) r dr)
 * overestimate the 2p3/2 -> 1s rate at high Z by ~20% because the photon
 * momentum k = omega/c is not small compared with the momentum of the
 * inner electron (k ~ 0.2 a.u. at Z = 92).  The fully retarded rate is
 * therefore computed from the four-component plane-wave matrix element
 *   M = <a | alpha . eps exp(-ik.r) | b>
 * expanded in partial waves j_L(kr) (L = 0..l1+l2), summed over all
 * magnetic quantum numbers and both photon polarisations:
 *   A = (2 k / (2 j1 + 1)) sum_{ma,mb,eps} |M|^2    (a.u.)
 * with k = omega/c the photon wavenumber.  No parity or multipole filter
 * is applied, so the result contains all multipoles allowed by angular
 * momentum (E1, M2, E3 ... for parity-changing pairs; M1, E2, M3 ... for
 * parity-conserving pairs).
 *
 * This reproduces the hydrogenic E1 rates of Gonzalez, Alexander &
 * Coldwell, J. Math. Chem. 61, 193 (2023), Table 1 for j = 1/2 rows
 * exactly and for j = 3/2 rows to ~1%, and the M1/E2 rates of Sapirstein,
 * Pachucki & Cheng (Phys. Rev. A 69, 022113 (2004)) and Knyazeva,
 * Lyashchenko, Zhang, Yu & Andreev (Phys. Rev. A 106, 012809 (2022)).
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the initial state
 * @param kappa1 angular quantum number of the initial state
 * @param n2     principal quantum number of the final state
 * @param kappa2 angular quantum number of the final state
 *
 * @return the fully retarded rate in s^-1, 0 if the transition is
 *         disallowed or has negative photon energy
 */

/**
 * One component of a spinor harmonic chi_{kappa,m}: (coeff, l, m).
 * (AngComp, spinor_harm, ang_int and sbesj are shared with the retarded RR
 * code in DiracRR.cxx and live in DiracAngle.h.)
 */

//////////////////////////////////////////////////////////////////////////
static double retarded_total(double z, int n1, int kappa1, int n2, int kappa2) {
  using hydroconst::alpha;
  using hydroconst::au_t_s;

  // Reject non-existent states (l >= n)
  if (!valid_kappa(n1, kappa1) || !valid_kappa(n2, kappa2)) return 0.0;

  int l1, l2;
  double j1, j2;
  lj_from_kappa(kappa1, l1, j1);
  lj_from_kappa(kappa2, l2, j2);

  // Delta j <= 2 with the E2 exclusions (0->0, 0->1 and 1->0).  No parity
  // filter: all multipoles allowed by angular momentum contribute.
  if (fabs(j1 - j2) > 2.01) return 0.0;
  if (j1 < 0.01 && j2 < 0.01) return 0.0;
  if (j1 < 0.01 && fabs(j2 - 1.0) < 0.01) return 0.0;
  if (j2 < 0.01 && fabs(j1 - 1.0) < 0.01) return 0.0;

  double c = 1.0 / alpha;

  // Photon wavenumber in atomic units (initial state 1, final state 2)
  double E1 = dirac_binding_energy(z, n1, kappa1);       // eV, initial
  double E2 = dirac_binding_energy(z, n2, kappa2);       // eV, final
  double omega = (E1 - E2) / (2.0 * hydroconst::Ryd_eV); // Hartree
  if (omega <= 0.0) return 0.0;
  double k = omega / c;

  // Log-spaced radial grid (Bohr radii), as in the reference code
  const int ng = 160000;
  const double r_min = 1e-10, r_max = 80.0;
  const double dlnr = log(r_max / r_min) / (ng - 1);

  // Largest angular momentum difference needed by the spinor-harmonic pairs
  auto lof = [&](int kp) {
    int l;
    double j;
    lj_from_kappa(kp, l, j);
    return l;
  };
  int Lmax = max(lof(kappa2) + lof(-kappa1), lof(-kappa2) + lof(kappa1));

  // Radial wavefunctions on the grid; P (large) and Q (small) with
  // int (P^2 + Q^2) dr = 1.
  vector<double> r(ng), Pa(ng), Qa(ng), Pb(ng), Qb(ng);
  for (int i = 0; i < ng; i++) r[i] = r_min * exp(dlnr * i);

  double za = z * alpha;
  auto radial_state = [&](int n, int kappa, vector<double> &P, vector<double> &Q) {
    int nr = n - abs(kappa);
    double gamma = sqrt(kappa * kappa - za * za);
    double Nprime = nr + gamma;
    double eps = 1.0 / sqrt(1.0 + za * za / (Nprime * Nprime));
    double lam = z / sqrt(Nprime * Nprime + za * za);
    double two_lam = 2.0 * lam;
    double C = dirac_normalisation(z, n, kappa);
    vector<double> acoef, bcoef;
    build_FG_coeffs(kappa, nr, gamma, eps, two_lam, z, acoef, bcoef);
    double gp = C * sqrt(1.0 + eps);
    double gm = C * sqrt(1.0 - eps);
    for (int i = 0; i < ng; i++) {
      double rp = pow(two_lam * r[i], gamma) * exp(-lam * r[i]);
      double F = 0.0, G = 0.0, rk = 1.0;
      for (int jj = 0; jj <= nr; jj++) {
        F += acoef[jj] * rk;
        G += bcoef[jj] * rk;
        rk *= r[i];
      }
      P[i] = gp * rp * F;
      Q[i] = gm * rp * G;
    }
  };

  radial_state(n2, kappa2, Pa, Qa); // final (a)
  radial_state(n1, kappa1, Pb, Qb); // initial (b)

  // Spherical Bessel functions j_L(kr) for L = 0..Lmax
  vector<vector<double> > J(Lmax + 1, vector<double>(ng));
  for (int i = 0; i < ng; i++) {
    double x = k * r[i];
    if (x < 1e-8) {
      J[0][i] = 1.0;
      for (int L = 1; L <= Lmax; L++) J[L][i] = 0.0;
      continue;
    }
    for (int L = 0; L <= Lmax; L++) J[L][i] = sbesj(L, x);
  }

  // Radial integrals with the small/large cross combinations
  //   Rp[L] = int Pa Qb j_L(kr) dr      (a large, b small)
  //   Rm[L] = int Qa Pb j_L(kr) dr      (a small, b large)
  vector<double> Rp(Lmax + 1, 0.0), Rm(Lmax + 1, 0.0);
  for (int L = 0; L <= Lmax; L++) {
    for (int i = 1; i < ng; i++) {
      double rp_prev = Pa[i - 1] * Qb[i - 1] * J[L][i - 1] * r[i - 1];
      double rp_cur = Pa[i] * Qb[i] * J[L][i] * r[i];
      double rm_prev = Qa[i - 1] * Pb[i - 1] * J[L][i - 1] * r[i - 1];
      double rm_cur = Qa[i] * Pb[i] * J[L][i] * r[i];
      Rp[L] += 0.5 * (rp_prev + rp_cur);
      Rm[L] += 0.5 * (rm_prev + rm_cur);
    }
    Rp[L] *= dlnr;
    Rm[L] *= dlnr;
  }

  // Sum over all magnetic quantum numbers of initial (b) and final (a)
  int ma_min = -(int)(2.0 * j2), ma_max = (int)(2.0 * j2);
  int mb_min = -(int)(2.0 * j1), mb_max = (int)(2.0 * j1);
  double tot = 0.0;
  for (int ms2a = ma_min; ms2a <= ma_max; ms2a += 2) {
    AngComp A0, A1, A2, A3;
    spinor_harm(kappa2, ms2a, A0, A1);
    spinor_harm(-kappa2, ms2a, A2, A3);
    for (int ms2b = mb_min; ms2b <= mb_max; ms2b += 2) {
      AngComp B0, B1, B2, B3;
      spinor_harm(kappa1, ms2b, B0, B1);
      spinor_harm(-kappa1, ms2b, B2, B3);

      // component pairs (large a x small b: +i, small a x large b: -i)
      //   pair 0: (A0, B3)  pair 1: (A1, B2)
      //   pair 2: (A2, B1)  pair 3: (A3, B0)
      // Mx = i  [(A0B3) + (A1B2) - (A2B1) - (A3B0)]
      // My =     [(A0B3) - (A1B2) - (A2B1) + (A3B0)]
      complex<double> S0(0.0, 0.0), S1(0.0, 0.0), S2(0.0, 0.0), S3(0.0, 0.0);
      {
        vector<complex<double> > ac = ang_int(A0.l, A0.m, B3.l, B3.m);
        for (size_t L = 0; L < ac.size(); L++)
          if (ac[L] != complex<double>(0.0, 0.0)) S0 += ac[L] * Rp[L];
      }
      {
        vector<complex<double> > ac = ang_int(A1.l, A1.m, B2.l, B2.m);
        for (size_t L = 0; L < ac.size(); L++)
          if (ac[L] != complex<double>(0.0, 0.0)) S1 += ac[L] * Rp[L];
      }
      {
        vector<complex<double> > ac = ang_int(A2.l, A2.m, B1.l, B1.m);
        for (size_t L = 0; L < ac.size(); L++)
          if (ac[L] != complex<double>(0.0, 0.0)) S2 += ac[L] * Rm[L];
      }
      {
        vector<complex<double> > ac = ang_int(A3.l, A3.m, B0.l, B0.m);
        for (size_t L = 0; L < ac.size(); L++)
          if (ac[L] != complex<double>(0.0, 0.0)) S3 += ac[L] * Rm[L];
      }

      complex<double> c03 = A0.coeff * B3.coeff * S0;
      complex<double> c12 = A1.coeff * B2.coeff * S1;
      complex<double> c21 = A2.coeff * B1.coeff * S2;
      complex<double> c30 = A3.coeff * B0.coeff * S3;
      complex<double> Mx = complex<double>(0.0, 1.0) * (c03 + c12 - c21 - c30);
      complex<double> My = c03 - c12 - c21 + c30;
      tot += norm(Mx) + norm(My);
    }
  }

  // A = 2k / (2j1+1) * sum |M|^2  in atomic units
  double A_au = 2.0 * k / (2.0 * j1 + 1.0) * tot;
  return A_au / au_t_s;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Dirac E1 transition rate in s^-1
 *
 * E1 selection: parity must change (l1 + l2 odd) and Delta j <= 1.
 * kappa2 = -kappa1 (e.g. s1/2 <-> p1/2) and kappa2 = kappa1 +/- 1
 * (e.g. p1/2 <-> d3/2) are both allowed.
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the initial state
 * @param kappa1 angular quantum number of the initial state
 * @param n2     principal quantum number of the final state
 * @param kappa2 angular quantum number of the final state
 *
 * @return the E1 rate in s^-1, 0 if the transition is disallowed
 */
double dirac_transrate(double z, int n1, int kappa1, int n2, int kappa2) {
  if (!valid_kappa(n1, kappa1) || !valid_kappa(n2, kappa2)) return 0.0;

  int l1, l2;
  double j1, j2;
  lj_from_kappa(kappa1, l1, j1);
  lj_from_kappa(kappa2, l2, j2);

  // E1 selection: parity must change (l1 + l2 odd) and Delta j <= 1.
  // kappa2 = -kappa1 (e.g. s1/2 <-> p1/2) and kappa2 = kappa1 +/- 1
  // (e.g. p1/2 <-> d3/2) are both allowed.
  if ((l1 + l2 + 1) % 2 != 0) return 0.0;
  if (fabs(j1 - j2) > 1.01) return 0.0;
  if (j1 < 0.01 && j2 < 0.01) return 0.0;

  return retarded_total(z, n1, kappa1, n2, kappa2);
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Dirac E2 transition rate in s^-1
 *
 * E2: parity conserving (l1 + l2 even) with Delta l = +/- 2 so that E2 is
 * the leading multipole, and Delta j <= 2 with the forbidden 0 -> 0, 0 -> 1
 * and 1 -> 0 cases.  The rate is the fully retarded plane-wave rate, which
 * for such pairs is dominated by E2 (higher multipoles are (alpha Z)^2
 * suppressed).  Pairs with Delta l = 0 are handled by dirac_m1_rate.
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the initial state
 * @param kappa1 angular quantum number of the initial state
 * @param n2     principal quantum number of the final state
 * @param kappa2 angular quantum number of the final state
 *
 * @return the E2 rate in s^-1, 0 if the transition is disallowed
 */
double dirac_e2_rate(double z, int n1, int kappa1, int n2, int kappa2) {
  if (!valid_kappa(n1, kappa1) || !valid_kappa(n2, kappa2)) return 0.0;

  int l1, l2;
  double j1, j2;
  lj_from_kappa(kappa1, l1, j1);
  lj_from_kappa(kappa2, l2, j2);

  // E2: parity conserving (l1 + l2 even) with Delta l = +/- 2 so that E2 is
  // the leading multipole, and Delta j <= 2 with the forbidden 0 -> 0, 0 -> 1
  // and 1 -> 0 cases.  The rate is the fully retarded plane-wave rate, which
  // for such pairs is dominated by E2 (higher multipoles are (alpha Z)^2
  // suppressed).  Pairs with Delta l = 0 are handled by dirac_m1_rate.
  if ((l1 + l2) % 2 != 0) return 0.0;
  if (abs(l1 - l2) != 2) return 0.0;
  if (fabs(j1 - j2) > 2.01) return 0.0;
  if (j1 < 0.01 && j2 < 0.01) return 0.0;
  if (j1 < 0.01 && fabs(j2 - 1.0) < 0.01) return 0.0;
  if (j2 < 0.01 && fabs(j1 - 1.0) < 0.01) return 0.0;

  return retarded_total(z, n1, kappa1, n2, kappa2);
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief fine-structure M1 (magnetic dipole) rate in s^-1
 *
 * Magnetic-dipole decay between the two fine-structure partners of the same
 * (n, l) doublet, j = l+1/2 -> j' = l-1/2 (Delta l = 0, Delta j = 1).  The
 * dipole operator is mu = mu_B (L + 2S); J connects only equal-j states, so
 * the transition is driven by the spin S alone.  The reduced element is
 *   |<j'||L+2S||j>|^2 / (2j+1) = l / (2l+1).
 * The fully retarded plane-wave calculation (retarded_total) loses this
 * tiny rate to catastrophic cancellation at the very small photon momentum
 * of the near-degenerate fine-structure splitting, so the standard
 * magnetic-dipole formula is used here (SI units):
 *   A = (4/3) omega^3 mu_B^2 / (hbar c^3) * l / (2l+1).
 * It reproduces the hydrogen 2p3/2 -> 2p1/2 rate A ~ 4.4e-6 s^-1
 * (tau ~ 2.6 days) using the Dirac fine-structure splitting of 10.97 GHz.
 *
 * This formula assumes the two fine-structure partners share essentially
 * the same (non-relativistic) radial wavefunction, which only holds for
 * small Z*alpha. Extrapolated to high Z it diverges as omega^3 with the
 * rapidly growing Dirac splitting and becomes unphysical -- e.g. at Z=82 it
 * exceeds the fully allowed 2p3/2->1s1/2 E1 rate by a factor of ~34, and at
 * Z=92 by a factor of ~113, which cannot be correct for a magnetic-dipole
 * "forbidden" line. It is therefore restricted to Z <= FSM1_ZMAX below,
 * comfortably inside the regime it was validated for; heavier ions fall
 * through to retarded_total (dirac_m1_rate), which gives a small, physically
 * sensible M1/E1 ratio once Z is large enough that its own near-degenerate
 * cancellation problem no longer applies.
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the upper state
 * @param kappa1 angular quantum number of the upper state
 * @param n2     principal quantum number of the lower state
 * @param kappa2 angular quantum number of the lower state (fine-structure
 *                partner of kappa1)
 *
 * @return the M1 fine-structure rate in s^-1, 0 if the splitting is not
 *         positive
 */
// Above this nuclear charge the near-degenerate approximation behind
// fine_structure_m1 is no longer valid (see docstring above); dirac_m1_rate
// switches to retarded_total instead. Chosen from direct comparison: with
// this diagnostic build, retarded_total is itself numerically unreliable
// (erratic, non-monotonic in Z) for Z <~ 4-5 but settles into a smooth,
// physically sensible (M1 << E1) trend from Z = 5 upward.
static const int FSM1_ZMAX = 4;

static double fine_structure_m1(double z, int n1, int kappa1, int n2,
                                int kappa2) {
  using hydroconst::Ryd_eV;
  using hydroconst::au_t_s;

  int l;
  double j;
  lj_from_kappa(kappa1, l, j); // common l (kappa2 is the fine-structure partner)

  double E_up = dirac_binding_energy(z, n1, kappa1);
  double E_lo = dirac_binding_energy(z, n2, kappa2);
  double omega_au = (E_up - E_lo) / (2.0 * Ryd_eV); // Hartree
  if (omega_au <= 0.0) return 0.0;

  double omega_s = omega_au / au_t_s; // 1/s

  const double muB = 9.2740100783e-24; // J/T (CODATA 2018)
  const double hbar = 1.054571817e-34; // J s (CODATA 2018)
  const double c = 299792458.0;        // m/s (CODATA 2018)

  double B2 = l / (2.0 * l + 1.0); // |<j'||L+2S||j>|^2 / (2j+1)
  return (4.0 / 3.0) * omega_s * omega_s * omega_s * muB * muB /
         (hbar * c * c * c) * B2;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Dirac M1 transition rate in s^-1
 *
 * M1: Delta l = 0 and Delta j = 0, +/- 1 (0 -> 0 forbidden).
 *
 * Fine-structure decay within the same (n, l) doublet uses the analytic
 * magnetic-dipole formula (see fine_structure_m1) only for Z <= FSM1_ZMAX,
 * where it is validated and retarded_total is not.  All other Delta l = 0
 * pairs -- including same-(n,l) doublets at higher Z, and pairs such as
 * 2s -> 1s driven by the small components of the wavefunctions -- use the
 * fully retarded plane-wave rate.
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the initial state
 * @param kappa1 angular quantum number of the initial state
 * @param n2     principal quantum number of the final state
 * @param kappa2 angular quantum number of the final state
 *
 * @return the M1 rate in s^-1, 0 if the transition is disallowed
 */
double dirac_m1_rate(double z, int n1, int kappa1, int n2, int kappa2) {
  // Reject non-existent states (l >= n)
  if (!valid_kappa(n1, kappa1) || !valid_kappa(n2, kappa2)) return 0.0;

  int l1, l2;
  double j1, j2;
  lj_from_kappa(kappa1, l1, j1);
  lj_from_kappa(kappa2, l2, j2);

  // M1: Delta l = 0 and Delta j = 0, +/- 1 (0 -> 0 forbidden).
  if (l1 != l2) return 0.0;
  if (fabs(j1 - j2) > 1.01) return 0.0;
  if (j1 < 0.01 && j2 < 0.01) return 0.0;

  // Fine-structure decay within the same (n, l) doublet: use the analytic
  // magnetic-dipole formula only where it is valid (see above); otherwise
  // fall through to the fully retarded plane-wave rate below.
  if (n1 == n2 && kappa1 != kappa2 && z <= FSM1_ZMAX)
    return fine_structure_m1(z, n1, kappa1, n2, kappa2);

  return retarded_total(z, n1, kappa1, n2, kappa2);
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief total transition rate including E1, E2, M1 in s^-1
 *
 * Returns 0 when the transition is forbidden by all multipole selection
 * rules.
 *
 * E1 requires a parity change (l1+l2 odd), whereas E2 requires
 * l1+l2 even and M1 requires l1==l2, so E1 is always mutually
 * exclusive with E2 and M1 for a given pair; summing all three
 * multipoles therefore never double-counts.
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the initial state
 * @param kappa1 angular quantum number of the initial state
 * @param n2     principal quantum number of the final state
 * @param kappa2 angular quantum number of the final state
 *
 * @return the total rate in s^-1, 0 if the transition is forbidden
 */
double dirac_total_transrate(double z, int n1, int kappa1,
                              int n2, int kappa2) {
  if (!valid_kappa(n1, kappa1) || !valid_kappa(n2, kappa2)) return 0.0;
  // E1 requires a parity change (l1+l2 odd), whereas E2 requires
  // l1+l2 even and M1 requires l1==l2, so E1 is always mutually
  // exclusive with E2 and M1 for a given pair; summing all three
  // multipoles therefore never double-counts.
  return dirac_transrate(z, n1, kappa1, n2, kappa2)
       + dirac_e2_rate(z, n1, kappa1, n2, kappa2)
       + dirac_m1_rate(z, n1, kappa1, n2, kappa2);
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief lifetime (s) including E1, E2, and M1 transitions
 *
 * @param z     nuclear charge
 * @param n     principal quantum number of the level
 * @param kappa relativistic angular quantum number of the level
 * @param nmax  maximum principal quantum number considered (unused,
 *              decays are summed up to n only)
 *
 * @return the lifetime in seconds, or 10 s if no allowed decay exists
 *         (metastable fallback)
 */
double dirac_lifetime(double z, int n, int kappa, int /*nmax*/) {
  double sum = 0.0;

  for (int np = 1; np <= n; np++) {
    for (int kappap = -np; kappap <= np; kappap++) {
      if (kappap == 0) continue;
      if (!valid_kappa(np, kappap)) continue;
      if (np == n && kappap == kappa) continue;
      double Eup = dirac_binding_energy(z, n, kappa);
      double Elo = dirac_binding_energy(z, np, kappap);
      if (Eup <= Elo) continue;
      sum += dirac_total_transrate(z, n, kappa, np, kappap);
    }
  }
  return sum > 0.0 ? 1.0 / sum : 10.0;
}




//////////////////////////////////////////////////////////////////////////
/**
 * @brief branching ratio of a Dirac E1 transition
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the initial state
 * @param kappa1 angular quantum number of the initial state
 * @param n2     principal quantum number of the final state
 * @param kappa2 angular quantum number of the final state
 * @param nmax   maximum principal quantum number used in the lifetime
 *
 * @return the branching ratio tau * A (dimensionless)
 */
double dirac_branch(double z, int n1, int kappa1, int n2, int kappa2,
                     int nmax) {
  double tau = dirac_lifetime(z, n1, kappa1, nmax);
  return tau * dirac_transrate(z, n1, kappa1, n2, kappa2);
}
