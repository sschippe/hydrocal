/**
    @file DiracTrans.cxx
    @brief Calculation of relativistic (Dirac) hydrogenic transition rates
    and lifetimes (counterpart of lifetime.cxx)

    States are specified by (n, kappa); kappa encodes l and j:
    kappa = -(l+1) for j = l+1/2,  kappa = +l for j = l-1/2.

    @par CREATION
    @author Stefan Schippers
    @date 2026

    @par VERSION
    @verbatim
    $Id: DiracTrans.cxx 2098 2026-08-06 13:21:05Z iamp $
// SPDX-License-Identifier: MIT
    @endverbatim

 */
#include "DiracTrans.h"
#include "DiracRate.h"
#include "buildinfo.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

// ---------------------------------------------------------------------------
// The three basic quantities (simple wrappers over the diracrate library)
// ---------------------------------------------------------------------------

double diractrans(double z, int n1, int kappa1, int n2, int kappa2) {
  return dirac_total_transrate(z, n1, kappa1, n2, kappa2);
}

double diraclife(double z, int n, int kappa) {
  double sum = 0.0;
  for (int np = 1; np <= n; np++) {
    for (int kappap = -np; kappap <= np; kappap++) {
      if (kappap == 0) continue;
      if (np == n && kappap == kappa) continue;
      if (dirac_binding_energy(z, np, kappap) >=
          dirac_binding_energy(z, n, kappa))
        continue;
      sum += diractrans(z, n, kappa, np, kappap);
    }
  }
  return sum > 0.0 ? 1.0 / sum : 10.0; // 10 s = metastable fallback
}

double diracbranch(double z, int n1, int kappa1, int n2, int kappa2) {
  return diraclife(z, n1, kappa1) * diractrans(z, n1, kappa1, n2, kappa2);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Format a state as "n,kappa" with its (l,j) attached, e.g. "(2,-2)[1,1.5]"
static void statelabel(char *buf, size_t len, int n, int kappa) {
  int l;
  double j;
  lj_from_kappa(kappa, l, j);
  snprintf(buf, len, "(%d,%+d)[%d,%.1f]", n, kappa, l, j);
}

static int state_valid(int n, int kappa) {
  int l;
  double j;
  lj_from_kappa(kappa, l, j);
  return (l < n) ? 1 : 0;
}

// ---------------------------------------------------------------------------
// Menu-option-2 style output functions
// ---------------------------------------------------------------------------

// All n -> n' transitions for fixed n1, n2 (upper n1, lower n2)
static void dirac_fixed_n1_n2(void) {
  string filename;
  double z, tr, lt;
  int n1, n2, kappa1, kappa2;
  printf("\n Give Z, n and n' for n -> n' transition : ");
  scanf("%lf %d %d", &z, &n1, &n2);
  printf("\n Give filename (*.dtr) for output ........: ");
  cin >> filename;

  ofstream fout;
  filename += ".dtr";
  fout.open(filename);

  fout << "####################################################################################" << endl;
  fout << "### Fully relativistic hydrogenic transition rates acccounting for full retardation"  << endl;
  fout << "###" << endl;
  fout << "###  hydrocal revision     : " << HYDROCAL_REVISION
       << endl; // defined in buildinfo.h
  fout << "###               filename : " << filename << endl;
  fout << "###" << endl;
  fout << "###         nuclear charge : " << z << endl;
  fout << "###                upper n : " << n1 << endl;
  fout << "###                lower n : " << n2 << endl;
  fout << "###" << endl;
  fout << "###        upper state  ->  lower state       ";
  fout << ":   life-        transition     branching          E1           E2           M1" << endl;
  fout << "###      n,kappa [l,j]  ->  n',kappa' [l',j'] ";
  fout << ":   time[s]      rate [1/s]     ratio            [1/s]        [1/s]        [1/s]" << endl;
  fout << "###-----------------------------------------------------------------------------------------------------------------------------" << endl;
  for (kappa1 = -n1; kappa1 <= n1; kappa1++) {
    if (kappa1 == 0) continue;
    if (!state_valid(n1, kappa1)) continue;
    lt = diraclife(z, n1, kappa1);
    for (kappa2 = -n2; kappa2 <= n2; kappa2++) {
      if (kappa2 == 0) continue;
      if (!state_valid(n2, kappa2)) continue;
      if (dirac_binding_energy(z, n2, kappa2) >=
          dirac_binding_energy(z, n1, kappa1))
        continue;
      double A_E1 = dirac_transrate(z, n1, kappa1, n2, kappa2);
      double A_E2 = dirac_e2_rate(z, n1, kappa1, n2, kappa2);
      double A_M1 = dirac_m1_rate(z, n1, kappa1, n2, kappa2);
      tr = A_E1 + A_E2 + A_M1;
      if (tr <= 0.0) continue;
      double br = lt * tr;
      char b1[32], b2[32];
      statelabel(b1, sizeof(b1), n1, kappa1);
      statelabel(b2, sizeof(b2), n2, kappa2);
      fout << setprecision(5);
      fout << setw(22) << b1 << " ";
      fout << setw(22) << b2 << " : ";
      fout << setw(12) << lt << " ";
      fout << setw(12) << tr << " ";
      fout << setw(13) << br << " ";
      fout << setw(12) << A_E1 << " ";
      fout << setw(12) << A_E2 << " ";
      fout << setw(12) << A_M1 << endl;
    }
  }
fout.close();
printf("\n Transition data written to file %s.\n", filename.c_str());
}

// Fixed lower state n2,kappa2 and upper l (kappa1); rates for n1 = n2+1..nmax
static void dirac_fixed_n2_k2_k1(void) {
  double z, tr;
  int n1, kappa1, kappa2, n2, nmax;
  printf("\n Give Z, n', kappa', nmax, kappa: ");
  scanf("%lf %d %d %d %d", &z, &n2, &kappa2, &nmax, &kappa1);
  if (!state_valid(n2, kappa2)) {
    printf("\n n',kappa' not a valid state.\n");
    return;
  }
  printf("\n  n, kappa ->  n', kappa' :   rate [1/s], rate*n^3 [1/s]\n");
  for (n1 = n2 + 1; n1 <= nmax; n1++) {
    if (!state_valid(n1, kappa1)) continue;
    tr = diractrans(z, n1, kappa1, n2, kappa2);
    if (tr <= 0.0) continue;
    char b1[32], b2[32];
    statelabel(b1, sizeof(b1), n1, kappa1);
    statelabel(b2, sizeof(b2), n2, kappa2);
    printf(" %14s -> %14s : %12.4g, %12.4g\n", b1, b2, tr,
           tr * n1 * n1 * n1);
  }
}

// All relativistic lifetimes up to nmax (output to *.tau file)
static void dirac_all_lifetimes(void) {
  double z;
  int nmax, choice;
  string filename;

  printf("\n Give nuclear charge Z and maximum n : ");
  scanf("%lf %d", &z, &nmax);
  printf("\n Which kind of output?");
  printf("\n 1: lifetimes in s");
  printf("\n 2: lifetimes relative to the maximum value per n");
  printf("\n 3: transition rates in 1/s");
  printf("\n Make a choice ......................: ");
  scanf("%d", &choice);
  printf("\n Give filename for output (*.tau) ...: ");
  cin >> filename;

  ofstream fout;
  filename += ".tau";
  fout.open(filename);

  for (int n = 1; n <= nmax; n++) {
    double taumax = 0.0;
    vector<double> tau;
    for (int kappa = -n; kappa <= n; kappa++) {
      if (kappa == 0) continue;
      if (!state_valid(n, kappa)) continue;
      double t = diraclife(z, n, kappa);
      tau.push_back(t);
      if (t > taumax) taumax = t;
    }
    for (size_t k = 0; k < tau.size(); k++) {
      switch (choice) {
      case 2:
        fout << " " << setw(12) << setprecision(6) << tau[k] / taumax;
        break;
      case 3:
        fout << " " << setw(12) << setprecision(6)
             << (tau[k] > 0 ? 1.0 / tau[k] : 999.0);
        break;
      default:
        fout << " " << setw(12) << setprecision(6) << tau[k];
        break;
      }
    }
    fout << "\n";
  }
  fout.close();
  printf("\n Written relativistic lifetimes/rates to %s.tau\n",
         filename.c_str());
}

// Lifetime of a specific n,kappa level
static void dirac_lifetime_nk(void) {
  double z;
  int n, kappa;
  printf("\n Give nuclear charge Z and n, kappa : ");
  scanf("%lf %d %d", &z, &n, &kappa);
  if (!state_valid(n, kappa)) {
    printf("\n n,kappa not a valid state.\n");
    return;
  }
  double tau = diraclife(z, n, kappa);
  printf("\n        Lifetime: %14.8g s", tau);
  printf("\n Transition rate: %14.8g /s", 1.0 / tau);
}

void testDiracLifetime(void) {
  int choice;
  printf("\n 1: Fixed n', kappa' and kappa");
  printf("\n 2: all n -> n' transitions");
  printf("\n 3: all relativistic lifetimes up to nmax");
  printf("\n 4: lifetime of a specific n,kappa level");
  printf("\n\n Make a choice ...........: ");
  scanf("%d", &choice);

  switch (choice) {
  case 2:
    dirac_fixed_n1_n2();
    break;
  case 3:
    dirac_all_lifetimes();
    break;
  case 4:
    dirac_lifetime_nk();
    break;
  default:
    dirac_fixed_n2_k2_k1();
    break;
  }
}
