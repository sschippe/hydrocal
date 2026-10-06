/**
 * @file DRmodelspe.cxx
 *
 * @brief Calculation of model DR spectra from parameterized atomic rates
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */
#include "CoolerFractions.h"
#include "clebsch.h"
#include "fele.h"
#include "hydroconst.h"
#include "lifetime.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
static const double clight = hydroconst::clight_cm_s;

void calcdr(void) {
  const double ryd = hydroconst::Ryd_eV; // eV
  const double mec2 = hydroconst::mec2_eV;       // eV
  const int maxn = 1000;                         // highest n+1

  char answer;
  double qdef[maxn], q, slim, arI, arII, aa0, ld1, ld2, ld3, ld4, factor,
      q2 = 1.0, q4 = 1.0;
  int model, hydrogenic = 0, field = 0, l=0, ncore=1, lcore=0, nstart = 1, lstart = 0,
             nmin=1, nmax=1, defmax=0;

  printf("\n Which model to use?");
  printf("\n       nl-dependent Auger rates (with DRF): 1");
  printf("\n   only n-dependent Auger rates   (no DRF): 2");
  printf("\n                            make a choice : ");
  scanf("%d", &model);

  printf("\n Give ion charge ..........................: ");
  scanf("%lf", &q);
  if (q < 1) {
    q = -q;
    q2 = q * q;
    q4 = q2 * q2;
  }

  printf("\n Give series limit in eV ..................: ");
  scanf("%lf", &slim);
  slim *= q2;
  if (model == 2) {
    printf("\n Give quantum defect ......................: ");
    scanf("%lf", &qdef[0]);
    defmax = 0;
  } else {
    l = 0;
    printf("\n Different quantum defects can be spefied");
    printf("\n for each l separately. From some l onwards");
    printf("\n all quantum defects will be the same. The");
    printf("\n input can then be terminated by giving -1.");
    printf("\n If the input for the l=0 quantum defect");
    printf("\n is -1, then all quantum defects will be 0.\n");
    do {
      printf("\n Give quantum defect for l =%2d (-1 quits) .: ", l);
      scanf("%lf", &qdef[l]);
    } while (qdef[l++] >= 0);
    defmax = l - 2;
  }
  if (defmax < 0) {
    qdef[0] = 0.0;
    defmax = 0;
  }

  double emin, emax, edelta, ktpar, ktperp;
  printf("\n Give energy range (min, max, delta in eV) : ");
  scanf("%lf %lf %lf", &emin, &emax, &edelta);
  emax *= q2; // q2 is 1 if on input q>0 or q^2 if on input q<0

  printf("\n Give ktpar and ktperp in meV .............: ");
  scanf("%lf %lf", &ktpar, &ktperp);
  ktperp *= 0.001;
  ktpar *= 0.001;
  int i, epts = 1 + int((emax - emin) / edelta);
  vector<double> energy(epts, 0.0);
  vector<double> esigma(epts, 0.0);
  vector<double> esigmaF(epts, 0.0);

new_model:

  printf("\n Give constant factor S0 in eV2 cm2 s .....: ");
  scanf("%lf", &factor);

  printf("\n Give core transition rate ArI in 1/s .....: ");
  scanf("%lf", &arI);
  arI *= q4; // q4 is 1 if on input q>0 or q^4 if on input q<0

  printf("\n Radiative rates of type II transitions are");
  printf("\n modelled as ArII/n^3");
  if (model == 1) {
    printf("\n If a negative value is given for ArII then");
    printf("\n hydrogenic rates will be used instead.");
    printf("\n  -ArII is then interpreted as the number");
    printf("\n of fully occupied subshells in the ionic");
    printf("\n core, e.g. core = 1s2 2s: ArII = -1.");
    printf("\n Give ArII = -0.5 if no subshell is fully");
    printf("\n occupied (bare or hydrogenlike core).\n");
    printf("\n Give ArII in 1/s (if <0 hydrogenic rates).: ");
  } else {
    printf("\n Give ArII in 1/s .........................: ");
  }

  scanf("%lf", &arII);
  arII *= q4; // q4 is 1 if on input q>0 or q^4 if on input q<0
  if ((arII < 0) && (model == 1)) {
    hydrogenic = 1;
    if (arII > -0.9) {
      ncore = 0;
      lcore = 0;
    } else {
      // -arII = 0.5*ncore*(ncore-1)+lcore
      // calculate ncore and lcore from -arII
      double shell = sqrt(-2 * arII + 0.25) - 0.5;
      ncore = int(shell + 0.1);
      lcore = int(-arII - ncore * (ncore + 1) / 2 - 1 + 0.1);
      if (lcore >= ncore - 1) {
        nstart = ncore + 1;
        lstart = 0;
      } else {
        nstart = ncore;
        lstart = lcore + 1;
      }
    }
  } // end (if ((arII<0) ...

  if (model == 1) {
    printf("\n Auger rates are modeled as");
    printf(" Aa(n,l) = Aa0*exp(d1*l+d2*l^2+d3*l^3+d4*l^4)/n^3");
    printf("\n Give Aa0 (in 1/s), d1, d2, d3, d4 .......: ");
    scanf("%lf %lf %lf %lf %lf", &aa0, &ld1, &ld2, &ld3, &ld4);
  } else {
    printf("\n Auger rates are modelled as Aa0/n^3 (up to lmax)");
    printf("\n Give Aa0 (in 1/s) ........................: ");
    scanf("%lf", &aa0);
    ld1 = 0.0;
    ld2 = 0.0;
    ld3 = 0.0;
    ld4 = 0.0;
  }

  printf("\n Give nmin ................................: ");
  scanf("%d", &nmin);

new_nmax:
  printf("\n Give nmax ................................: ");
  scanf("%d", &nmax);

  int nAa = 0;
  if ((model == 1) && (defmax == 0)) {
    printf("\n Include electric field mixing ? (y/n) ....: ");
    scanf(" %c", &answer);
    field = (answer == 'y');
  }
  if (field) {
    printf("\n Write field mixed Auger rates Aa(n,k,m) to file ?");
    printf("\n Give one n value larger than 0 ..........:  ");
    scanf("%d", &nAa);
  }

  int fracdim = nmax * (nmax + 1) / 2;
  vector<double> fraction(fracdim, 1.0);
  string header;
  printf("\n Read surviving fractions from file ? (y/n): ");
  scanf(" %c", &answer);
  if (answer == 'y')
    readfraction(fraction, header, nmax);

  ofstream fout, fAnkm;
  string filename, fn, fnAnkm;
  printf("\n Give filename for output (*.spe) .........: ");
  cin >> fn;
  filename = fn + ".str";
  fout.open(filename);

  fout << "kTpar                : " << setw(10) << fixed << setprecision(5)
       << ktpar * 1000 << " meV\n";
  fout << "kTperp               : " << setw(10) << setprecision(5)
       << ktperp * 1000 << " meV\n";
  fout << "ion charge state     : " << setw(10) << setprecision(5) << q << "\n";
  fout << "core radiative rate  : " << setw(10) << setprecision(4)
       << defaultfloat << arI << " /s \n";
  if (hydrogenic) {
    fout << "type II rad. rates   : hydrogenic\n";
    fout << "ion core filled up to: n = " << setw(2) << ncore
         << ", l = " << setw(2) << lcore << "\n";
    fout << "highest final shell  : n = " << setw(2) << nmin - 1 << "\n";
  } else {
    fout << "type II rad. rate    : " << setw(10) << setprecision(4) << arII
         << " /s \n";
  }
  fout << "Aa0                  : " << setw(10) << setprecision(4) << aa0
       << " /s \n";
  if (model == 1) {
    fout << "l-decay coefficient 1: " << setw(10) << setprecision(4) << ld1
         << "\n";
    fout << "l-decay coefficient 2: " << setw(10) << setprecision(4) << ld2
         << "\n";
    fout << "l-decay coefficient 3: " << setw(10) << setprecision(4) << ld3
         << "\n";
    fout << "l-decay coefficient 4: " << setw(10) << setprecision(4) << ld4
         << "\n";
  }
  fout << "factor               : " << setw(10) << setprecision(4) << factor
       << " eV2 cm2 s\n";
  fout << "series limit         : " << setw(10) << setprecision(5) << fixed
       << slim << " eV\n";
  for (l = 0; l < defmax; l++)
    fout << "q-defect for l  = " << setw(2) << l << " : " << setw(10)
         << setprecision(5) << fixed << qdef[l] << "\n";
  fout << "q-defect for l >= " << setw(2) << defmax << " : " << setw(10)
       << setprecision(5) << fixed << qdef[defmax] << "\n";
  fout << "\n\n In the follwing listed are:";
  fout << "\n       n: main quantum number";
  fout << "\n   E(eV): resonance energy in eV";
  fout << "\n      sn: DR resonance strength in eVcm^2";
  fout << "\n  Sum sn: sum of sn up to actual n";
  fout << "\n     nDR: number of states participating in DR";
  if (field) {
    fout << "\n     sFn: DR resonance strength at maximum field mixing";
    fout << "\n Sum sFn: sum of sFn up to actual n";
    fout << "\n    nDRF: number of mixed states participating in DR";
    fout << "\n    enh1: nDRF/nDR";
    fout << "\n    enh2: sn/sFn";
  }
  fout << "\n";

  fout << "\n    n      E(eV)         sn     Sum sn    nDR";
  if (field) {
    fout << "        sFn    Sum sFn   nDRF   enh1   enh2\n";
  } else {
    fout << "\n";
  }
  fout.close();

  // l-dependent part of the autoionization rates

  vector<double> lexpdecay(nmax, 0.0);

  int ldecayflag = 1;
  double ldecay = 0.0, ldecayold;
  for (l = 0; l < nmax; l++) {
    if (l > defmax)
      qdef[l] = qdef[defmax];
    ldecayold = ldecay;
    ldecay = (((ld4 * l + ld3) * l + ld2) * l + ld1) * l;
    if ((ldecayflag) && ((ldecay < -30) || ((l > 10) && (ldecay > ldecayold))))
      ldecayflag = 0;
    lexpdecay[l] = (ldecayflag) ? exp(ldecay) : -1.0;
  }

  if (nAa) { // ouput of Stark mixed autonionization rates for n = nAa
    fnAnkm = fn + ".Ankm";
    fAnkm.open(fnAnkm);

    double aa, j1 = 0.5 * (nAa - 1);
    for (int m = 0; m < nAa; m++) {
      int k;
      for (k = -nAa + 1; k < -nAa + m + 1; k++)
        fAnkm << " --";
      int kk = 1;
      for (k = -nAa + m + 1; k <= nAa - m - 1;
           k++) { // the quantum number k is counted in steps of 2
        if (kk) {
          aa = 0;
          for (l = m; l < nAa; l++) {
            if (lexpdecay[l] <= 0.0)
              break;
            double m1 = 0.5 * (m - k);
            double m2 = 0.5 * (m + k);
            double cg = CG(j1, j1, m1, m2, l, m);
            aa += cg * cg * lexpdecay[l];
          }
          aa *= aa0 / pow(double(nAa), 3);
          fAnkm << " " << setw(10) << setprecision(4) << aa;
          kk = 0;
        } else { // every other k does not occur
          fAnkm << " --";
          kk = 1;
        }
      } // end for (k...)
      for (k = nAa - m; k < nAa; k++)
        fAnkm << " --";
      fAnkm << "\n";
    } // end for (m...)
    fAnkm.close();
  } // end if(nAa)

  // initialize spectra
  for (i = 0; i < epts; i++) {
    energy[i] = emin + i * edelta;
    esigma[i] = 0.0;
    esigmaF[i] = 0.0;
  }

  // calculate spectra
  double aa, ar = 0.0, e = 0.0, es, lsum = 0.0, lsumF = 0.0,
             lsum2 = 0.0, nsum = 0.0, nsumF = 0.0;

  int arIIflag = hydrogenic;
  vector<double> arIIhydro;
  for (int n = nmin; n <= nmax; n++) {
    arIIhydro.resize(n);
    printf(" n = %4d", n);
    double n3 = n * n * n;
    double arIImax = 0; // maximum hydrogenic rate with one n manifold
    if (hydrogenic) {
      for (l = 0; l < n; l++) {
        arIIhydro[l] = 0.0;
        if (arIIflag) {
          for (int n2 = nstart; n2 < nmin; n2++) {
            for (int l2 = l - 1; l2 <= l + 1; l2 += 2) {
              if ((n2 == nstart) && (l2 < lstart))
                break;
              if ((l2 >= 0) && (l2 < n2))
                arIIhydro[l] += hydrotrans(q, n, l, n2, l2);
            }
          }
          if (arIIhydro[l] > arIImax)
            arIImax = arIIhydro[l];
        } // end if (arIIflag)
      } // end for (l...)
      if (arIImax < 0.01 * arI)
        arIIflag = 0;
    } else { // not hydrogenic
      ar = arI + arII / n3;
    }
    if (arIIflag)
      printf("   arIImax = %12.6g /s", arIImax);
    printf("\n");

    int countF = 0;
    int count0 = 0;
    if (model == 1) {
      lsum = 0;
      lsumF = 0;
      lsum2 = 0;
      for (l = 0; l < n; l++) {
        e = slim - ryd * q * q / (n - qdef[l]) / (n - qdef[l]);
        if (e < 0)
          continue;

        if (field) { // attention: l has the meaning of m now
          double j1 = 0.5 * (n - 1);
          for (int k = n - l - 1; k >= 0; k -= 2) {
            aa = 0;
            double arIIF = 0;
            for (int j = l; j < n;
                 j++) { // attention: j has the meaning of l now
              if (lexpdecay[j] <= 0.0)
                break;
              double l1 = 0.5 * (l - k);
              double l2 = 0.5 * (l + k);
              double cg = CG(j1, j1, l1, l2, j, l);
              aa += cg * cg * lexpdecay[j];
              if (arIIflag)
                arIIF += cg * cg * arIIhydro[j];
            }
            ar = arI + arIIF;
            aa *= aa0 / n3;
            int mult1 = l == 0 ? 1 : 2;
            int mult2 = k == 0 ? 1 : 2;
            lsumF += mult1 * mult2 * factor * ar * aa / (ar + aa);
            if (aa > ar)
              countF += mult1 * mult2;
          } // end for (k...)
        } // end if (field)
        if (lexpdecay[l] <= 0.0)
          break;
        aa = aa0 * lexpdecay[l] / n3;
        if (hydrogenic) {
          ar = arI + arIIhydro[l];
        } else {
          ar = arI + arII / n3;
        }
        if (aa > ar)
          count0 += 4 * l + 2;
        es = factor * (2 * l + 1) * ar * aa / (ar + aa) *
             fraction[(n - 1) * n / 2 + l];
        if (l <
            defmax) { // quantum defects are mutually different for l <= defmax
          for (i = 0; i < epts; i++) {
            esigma[i] += es * fecool(energy[i], e, ktpar, ktperp);
          }
          lsum2 += es;
        } // end if (l<defmax)
        lsum += es;
      } // end for(l...)
      for (i = 0; i < epts; i++) {
        double fe = fecool(energy[i], e, ktpar, ktperp);
        esigma[i] += (lsum - lsum2) * fe;
        if (field)
          esigmaF[i] += lsumF * fe;
      }
    } // end (model==1)
    else { // (model==2)
      e = slim - ryd * q * q / (n - qdef[0]) / (n - qdef[0]);
      aa = aa0 / n3;
      es = factor * aa * ar / (aa + ar);
      for (i = 0; i < epts; i++) {
        esigma[i] += es * fecool(energy[i], e, ktpar, ktperp);
      }
      lsum = es;
    } // end (model==2)
    nsum += lsum;
    nsumF += lsumF;
    if (e > 0) {
      fout.open(filename, ios::app);
      fout << setw(5) << n << " " << setw(10) << setprecision(4) << fixed << e
           << " " << setw(10) << setprecision(4) << defaultfloat << lsum << " "
           << setw(10) << setprecision(4) << nsum << " " << setw(6) << count0;
      if (field) {
        double enh1 = count0 > 0 ? double(countF) / (double(count0)) : 0.0;
        double enh2 = lsum > 0 ? lsumF / lsum : 0.0;
        fout << " " << setw(10) << setprecision(4) << lsumF << " " << setw(10)
             << setprecision(4) << nsumF << " " << setw(6) << countF << " "
             << setw(6) << setprecision(2) << fixed << enh1 << " " << setw(6)
             << setprecision(2) << fixed << enh2;
      }
      fout << "\n";
      fout.close();
    }
  } // end for (n...)
  cout << "\n List of line strengths written to file " << filename << "\n";

  filename = fn + ".spe";
  fout.open(filename);
  for (i = 0; i < epts; i++) {
    if (energy[i] == 0)
      continue;
    double alphafac = sqrt(2.0 / energy[i] / mec2) * clight;
    double alpha = esigma[i] > 1e-99 ? alphafac * esigma[i] : 0.0;
    fout << setw(12) << setprecision(6) << energy[i] << " " << setw(12)
         << setprecision(6) << alpha;
    if (field) {
      double alphaF = esigma[i] > 1e-99 ? alphafac * esigmaF[i] : 0.0;
      double enh = alpha > 0 ? alphaF / alpha : 0.0;
      fout << " " << setw(12) << setprecision(6) << alphaF << " " << setw(6)
           << setprecision(2) << fixed << enh;
    }
    fout << "\n";
  }
  fout.close();
  cout << "\n DR rate coefficient written to file " << filename << "\n";

  printf("\n New nmax ? (y/n) .........................: ");
  scanf(" %c", &answer);
  if (answer == 'y')
    goto new_nmax;

  printf("\n New model parameters ? (y/n) .............: ");
  scanf(" %c", &answer);
  if (answer == 'y')
    goto new_model;
}

double HFAugerfactor(double I, double j1core, double l1Rydberg, double F1,
                     double /*j2core*/, double F2) {
  // see Eq.9 of Pindzola et al., Phys. Rev. A 45 (1992) R7659
  double sum = 0.0;
  for (double j1 = fabs(l1Rydberg - 0.5); j1 <= l1Rydberg + 0.6; j1 += 1.0)
    for (double J = fabs(j1 - F1); J <= F1 + j1; J += 1.0)
      for (double j2 = fabs(F2 - J); j2 <= F2 + J + 0.1; j2 += 1.0) {
        //	printf("%4.1f %4.1f %4.1f\n",j1,J,j2);
        double ninej = NineJ(I, j1core, F1, j1core, l1Rydberg, j1, F2, j2, J);
        sum += (2 * j1 + 1) * (2 * j2 + 1) * ninej * ninej;
      }
  return 4 * sum * (2 * F1 + 1) * (2 * F2 + 1) * pow(2 * l1Rydberg + 1, -2);
}

void HyperfineDR(void) {
  double I, j1core, l1Rydberg, j2core;

  printf("\n Give I, j1core, l1Rydberg, j2core: ");
  scanf(" %lf %lf %lf %lf", &I, &j1core, &l1Rydberg, &j2core);

  double testsum = 0;
  for (double F1 = fabs(I - j1core); F1 <= I + j1core + 0.1; F1 += 1.0)
    for (double F2 = fabs(I - j2core); F2 <= I + j2core + 0.1; F2 += 1.0) {
      testsum += HFAugerfactor(I, j1core, l1Rydberg, F1, j2core, F2);
    }

  printf("\n testsum = %10.4g\n", testsum);
}
