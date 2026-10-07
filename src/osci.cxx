/**
 *   @file osci.cxx
 *
 *   @brief Calculation of hydrogenic bound-bound and bound-free oscillator strengths
 *
 *   @par CREATION
 *   @author Stefan Schippers
 *   @date 1997, 1999, 2023
 *
 *   @par VERSION
 *   @verbatim
// SPDX-License-Identifier: MIT
    @endverbatim
 *
 *
 *  Hydrogenic dipole oscillator strengths see
 *
 *  - Kazem Omidvar and Patricia T. Guimares,
 *        The Astrophysical Journal Suplement Series 73 (1990) 555-602,
 *  - H.A. Bethe and E.E. Salpeter, Quantum Mechanics of One and Two
 *         Electron Systems n Handbuch der Physik Vol XXXV (Springer, 1957).
 *
 * Calculations of radial hydrogenic matrix elements use recursion formulae
 * given by L. Infield and T.E.Hull, Rev. Mod. Phys 23 (1951) 21.
 *
*/
#include "hydromath.h"
#include "clebsch.h"
#include "osci.h"
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <vector>
#include "stdin_guard.h"

using namespace std;

//------------------------------------------------------------------------------------------------------
/**
 * @brief recursive calculation of bound-bound radial matrix elements for all l
 *
 * @param n1 principal quantum number of upper level
 * @param n2 principal quantum number of lower level
 * @param Rbbp on exit, array of matrix elements for l -> l+1
 * @param Rbbm on exit, array of matrix elements for l -> l-1
 */
void CalcRbb(int n1, int n2, vector<double> &Rbbp, vector<double> &Rbbm) {
  // start value
  int nm = n1 - n2;
  int np = n1 + n2;
  int n22 = n2 * 2;
  double nn2 = 2 * sqrt(double(n1 * n2));
  double xm = nm / nn2;
  double xp = nn2 / np;
  double fac = n1 * n2 * xp / double(np * np);

  if (nm == 1)
    fac *= nn2;

  for (int n = 0; n <= nm; n++) {
    if ((n >= 2) && (n < nm))
      fac *= xm / sqrt(double(n));
    int nii = n + n22;
    if ((nii >= n22) && (nii <= np))
      fac *= xp * sqrt(double(nii));
  }
  for (int n = 2; n < n22; n++)
    fac *= xp;
  Rbbm[n2] = fac;
  Rbbp[n2] = 0.0;

  // recursion
  double a, ap, am;
  for (int l = n2 - 1; l >= 0; l--) {
    int l1 = l + 1;
    a = 2 * sqrt(double((n2 + l) * (n2 - l))) / double(n2);
    am = (2 * l + 1) * sqrt(double((n1 + l1) * (n1 - l1))) / double(n1 * l1);
    ap = sqrt(double((n2 + l1) * (n2 - l1))) / double(n2 * l1);
    Rbbm[l] = (am * Rbbm[l1] + ap * Rbbp[l1]) / a;
    a = 2 * sqrt(double((n1 + l) * (n1 - l))) / double(n1);
    am = sqrt(double((n1 + l1) * (n1 - l1))) / double(n1 * l1);
    ap = (2 * l + 1) * sqrt(double((n2 + l1) * (n2 - l1))) / double(n2 * l1);
    Rbbp[l] = (am * Rbbm[l1] + ap * Rbbp[l1]) / a;
  }
}

//------------------------------------------------------------------------------------------------------
/**
 * @brief recursive calculation of bound-continuum radial matrix elements for
 * all l
 *
 * @param k wave number of the continuum electron
 * @param n2 principal quantum number of the bound electron
 * @param Rbcp on exit, array of matrix elements for l -> l+1
 * @param Rbcm on exit, array of matrix elements for l -> l-1
 */
void CalcRbc(double k, int n2, vector<double> &Rbcp, vector<double> &Rbcm) {
  // start values
  int n;
  double ksqr = k * k;
  double n2sqr = n2 * n2;
  double n2d = n2;
  double n22 = 2 * n2d;
  double x = 4.0 * n2d * k / (ksqr + n2sqr);
  double fac = sqrt(n2d * 0.5) * 0.25 * x * x * ksqr;
  for (n = 1; n <= n2; n++)
    fac *= x * sqrt((ksqr + n * n) / double(n * (n22 - n)));
  fac *= exp(-2.0 * k * atan(n2d / k));
  if (k < 10.0)
    fac /= sqrt(1.0 - exp(-8 * atan(1.0) * k));
  Rbcm[n2] = fac;
  Rbcp[n2] = 0.0;

  // recursion
  double a, ap, am;
  for (int l = n2 - 1; l >= 0; l--) {
    int l1 = l + 1;
    a = 2 * sqrt(double((n2 + l) * (n2 - l))) / double(n2);
    am = (2 * l + 1) * sqrt(k * k + l1 * l1) / double(k * l1);
    ap = sqrt(double((n2 + l1) * (n2 - l1))) / double(n2 * l1);
    Rbcm[l] = (am * Rbcm[l1] + ap * Rbcp[l1]) / a;
    a = 2 * sqrt(k * k + l * l) / double(k);
    am = sqrt(k * k + l1 * l1) / double(k * l1);
    ap = (2 * l + 1) * sqrt(double((n2 + l1) * (n2 - l1))) / double(n2 * l1);
    Rbcp[l] = (am * Rbcm[l1] + ap * Rbcp[l1]) / a;
  }
}

//------------------------------------------------------------------------------------------------------
/**
 * @brief hydrogenic bound-bound oscillator strength
 *
 * @param n1 principal quantum number of initial level
 * @param l1 orbital angular momentum quantum number of initial level
 * @param n2 principal quantum number of final level
 * @param l2 orbital angular momentum quantum number of final level
 *
 * @return bound-bound oscillator strength
 *
 * Since the radial matrix elements are calculated iteratively
 * all bound-bound matrix elements for given n1 and n2 are stored internally.
 * In order to avoid unnecessary calculations, in any summation of
 * oscillator strengths, sum over l before summing over n.
 */
double fosciBB(int n1, int l1, int n2, int l2) {
  static int n1_calculated = 0;
  static int n2_calculated = 0;
  static vector<double>
      Rbbp; // stores <n2, l1+1 | r | n1, l1> for fixed n1,n2, and all l1
  static vector<double>
      Rbbm; // stores <n2, l1-1 | r | n1, l1> for fixed n1 n2, and all l1

  if (n1 == n2)
    return 0.0;

  // radial matrix elements are only newly calculated when n1 or n2 have changed
  // from before
  if ((n1_calculated != n1) || (n2_calculated != n2)) {
    int ndim = n1 > n2 ? n1 + 1 : n2 + 1;
    Rbbp.resize(ndim, 0.0);
    Rbbm.resize(ndim, 0.0);
    if (n1 > n2) {
      CalcRbb(n1, n2, Rbbp, Rbbm);
    } else {
      CalcRbb(n2, n1, Rbbp, Rbbm);
    }
    n1_calculated = n1;
    n2_calculated = n2;
    // for (int nn=0; nn<ndim; nn++) printf("%g %g\n",Rbbm[nn],Rbbp[nn]);
  }

  double n1sqr = n1 * n1;
  double n2sqr = n2 * n2;
  double ll = 3.0 * (2.0 * l1 + 1.0);
  if (l2 == (l1 - 1)) {
    double r = n1 > n2 ? Rbbm[l1] : Rbbp[l2 + 1];
    return (1.0 / n1sqr - 1.0 / n2sqr) * l1 / ll * r * r;
  } else if (l2 == (l1 + 1)) {
    double r = n1 > n2 ? Rbbp[l1 + 1] : Rbbm[l2];
    return (1.0 / n1sqr - 1.0 / n2sqr) * l2 / ll * r * r;
  } else
    return 0;
}

//------------------------------------------------------------------------------------------------------
/**
 * @brief hydrogenic bound-free oscillator strength
 *
 * @param e continuum energy
 * @param n principal quantum number of the bound electron
 * @param l orbital angular momentum quantum number of the bound electron
 *
 * @return bound-free oscillator strength
 *
 * Since the radial matrix elements are calculated iteratively
 * all bound-bound matrix elements for given e and n are stored internally.
 * In order to avoid unnecessary calculations, in any summation of
 * oscillator strengths, sum over l before summing over n.
 */
double fosciBC(double e, int n, int l) {
  static int n_calculated = 0;
  static double k_calculated = -1.0;
  static vector<double> Rbcp, Rbcm;

  if (e <= 0.0)
    return 0.0;
  double k = 1.0 / sqrt(e);

  // radial matrix elements are only newly calculated when n or k have changed
  // from before
  if ((n_calculated != n) || (k_calculated != k)) {
    Rbcp.resize(n + 1, 0.0);
    Rbcm.resize(n + 1, 0.0);
    CalcRbc(k, n, Rbcp, Rbcm);
    k_calculated = k;
    n_calculated = n;
  }
  double fac = (e + 1.0 / n / n) / 3 / (2 * l + 1);
  double rm = Rbcm[l + 1];
  double rp = l > 0 ? Rbcp[l] : 0;

  return fac * ((l + 1) * rm * rm + l * rp * rp);
}


//------------------------------------------------------------------------------------------------------
/**
 * @brief bound-free matrix element squared for l->l+1 and zero energy
 *
 * @param n principal quantum number of bound electron
 * @param l orbital angular momentum quantum number of bound electron
 *
 * @return bound-free matrix element squared
 */
double rpsqr0(int n, int l) {
  double n4 = n * 4.0;
  double prod = pow(n4, 3);
  double ll = 0;
  for (int m = -l; m <= l; m++) {
    ll++;
    prod *= (n + m) * n4 / ll / ll;
  }
  double hyp1 = hypconfl(-n + l + 1, 2 * l + 2, n4);
  double hyp2 = hypconfl(-n + l + 2, 2 * l + 3, n4);
  return prod * pow(hyp1 + (n - l - 1.0) / (l + 1.0) * hyp2, 2);
}

//------------------------------------------------------------------------------------------------------
/**
 * @brief bound-free matrix element squared for l->l-1 and zero energy
 *
 * @param n principal quantum number of bound electron
 * @param l orbital angular momentum quantum number of bound electron
 *
 * @return bound-free matrix element squared
 */
double rmsqr0(int n, int l) {
  double n4 = n * 4.0;
  double prod = (n - l) * (n + l) * n4;
  double ll = 0;
  for (int m = -l + 1; m <= l - 1; m++) {
    ll++;
    prod *= (n + m) * n4 / ll / ll;
  }
  double hyp1 = hypconfl(-n + l + 1, 2 * l, n4);
  double hyp2 = hypconfl(-n + l - 1, 2 * l, n4);
  return prod * pow(hyp1 - hyp2, 2);
}

//------------------------------------------------------------------------------------------------------
/**
 * @brief hydrogenic bound-free oscillator strength at threshold
 *
 * @param n principal quantum number of bound electron
 * @param l orbital angular momentum quantum number of bound electron
 *
 * @return bound-free oscillator strength at threshold
 */
double fosciBC(int n, int l) {
  double r = l == 0 ? rpsqr0(n, l) : l * rmsqr0(n, l) + (l + 1) * rpsqr0(n, l);
  return r * exp(-4.0 * n) / (2 * l + 1) / 6.0;
}

//------------------------------------------------------------------------------------------------------
/**
 * @brief hydrogenic bound-free oscillator strength at threshold
 *
 * @param n principal quantum number of bound electron
 * @param l orbital angular momentum quantum number of bound electron
 * @param m magentic quantum number
 *
 * @return bound-free oscillator strength at threshold
 */
double fosciBC(int n, int l, int m) {
  double r;
  if (l == 0) {
    r = rpsqr0(n, 0);
  } else {
    double splus = 0.0, sminus = 0.0;
    for (int dm = -1; dm <= 1; dm++) {
      splus += Strength(l, l + 1, m, dm);
      sminus += Strength(l, l - 1, m, dm);
    }
    r = l * rmsqr0(n, l) * splus + (l + 1) * rpsqr0(n, l) * sminus;
  }
  return r * exp(-4.0 * n) / (2 * l + 1) / 6.0;
}

//------------------------------------------------------------------------------------------------------
/**
 * @brief test subroutine for bound-bound oscillator strength
 */
void testOsciBB(void) {
  int n1, l1, n2, l2;

  printf(" Give n,l and n',l' for n,l -> n',l' transition : ");
  scanf("%d %d %d %d", &n1, &l1, &n2, &l2);
  printf("\n\n fBB = %12.4G \n\n", fosciBB(n1, l1, n2, l2));
}

//------------------------------------------------------------------------------------------------------
/**
 * @brief test subroutine for bound-free oscillator strength
 */
void testOsciBC(void) {
  double x, f;
  int n, l;

  printf(" Give x (n*n*Rydberg), n, l : ");
  scanf("%lf %d %d", &x, &n, &l);
  f = x == 0 ? fosciBC(n, l) : fosciBC(x / n / n, n, l);
  printf("\n\n nfBC = %g /Ryd \n\n", f / n);
}
