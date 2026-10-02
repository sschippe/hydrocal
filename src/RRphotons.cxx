/**
 * @file RRphotons.cxx
 *
 * @brief Hydrogenic calculation of RR photon spectra including cascades
 *        Optionally uses relativistic (Dirac) transition rates for
 *        low principal quantum numbers.
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: RRphotons.cxx 2039 2026-07-20 07:57:32Z iamp $
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */
#include "RRphotons.h"
#include "RRDRxsec.h"
#include "DiracRate.h"
#include "DiracRR.h"
#include "hydroconst.h"
#include "hydromath.h"
#include "radrate.h"
#include "buildinfo.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
using hydroconst::pi;
//////////////////////////////////////////////////////////////////////////
/**
 * @brief Dirac energy of a n,l-level; weighted average over j-components
 *
 * @return energy in eV
 */
static double diracEnergy(double z, int n, int l) {
  double k1 = (double)l;
  double k2 = (double)(l + 1);
  double E1 = 0.0;
  double E2 = dirac_binding_energy(z, n, k2);  
  if (l == 0) {
    return E2; // only j = l+1/2 exists
  } else {
    E1 = dirac_binding_energy(z, n, k1);
    double Emean = 2.0*l*E1 + 2.0*(l + 1.0)*E2;
    Emean /= (double)(4.0*l+2.0);
    return Emean;
  }
}

// =========================================================================
// Relativistic / mixed cascade
// =========================================================================

struct StateInfo {
  int n, l;
  double j;      // 0 for unresolved (nonrelativistic)
  int kappa;     // 0 for unresolved
  bool is_rel;   // true = Dirac j-resolved
};

/////////////////////////////////////////////////////////////////////////////////
/**
 * @brief computes the number of states in mixed realtivistic and nonrelativistic cascade
 *
 * @param nmax maximum principal quantum number to be considered in the cascade computation
 * @param nrel n-value up to (inclusive) which levels and rates are treated relativistically
 *
 * @return number of levels in the cascade computation
 */
static int nstates_total(int nmax, int nrel) {
  if (nrel < 1) return nmax * (nmax + 1) / 2;
  int nnonrel = nmax * (nmax + 1) / 2 - nrel * (nrel + 1) / 2;
  return nrel*nrel + nnonrel;
}

/////////////////////////////////////////////////////////////////////////////////
/**
 * @brief gathers the relevant information about all levels to be considered
 *
 * @param states on exit, vector containg all required information
 * @param nmax maximum principal quantum number to be considered in the cascade computation
 * @param nrel n-value up to (inclusive) which levels and rates are treated relativistically
 */
static void build_state_info(vector<StateInfo> &states, int nmax, int nrel) {
  states.clear();
  // Relativistic states (n <= nrel) ordered by increasing |kappa|
  // so that state index increases with energy within each n.
  for (int n = 1; n <= nrel; n++) {
    for (int absk = 1; absk <= n; absk++) {
      // κ = -absk (l = absk-1, j = absk-1/2) — lower l, lower energy
      if (absk - 1 < n) {
        StateInfo s;
        s.n = n; s.is_rel = true; s.kappa = -absk;
        lj_from_kappa(s.kappa, s.l, s.j);
        states.push_back(s);
      }
      // κ = +absk (l = absk, j = absk-1/2)
      if (absk < n) {
        StateInfo s;
        s.n = n; s.is_rel = true; s.kappa = absk;
        lj_from_kappa(s.kappa, s.l, s.j);
        states.push_back(s);
      }
    }
  }
  // Nonrelativistic states (n > nrel)
  for (int n = nrel + 1; n <= nmax; n++) {
    for (int l = 0; l < n; l++) {
      StateInfo s;
      s.n = n;
      s.l = l;
      s.j = 0;
      s.kappa = 0;
      s.is_rel = false;
      states.push_back(s);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
/**
 * @brief Computation of the photon spectrum resulting from RR of a bare ion
 *
 * @param hydro object containing procomputed nonrelativistic hydrogenic transition rates
 * @param nmax maximum principal qauntum number considered in the cascade
 * @param z nuclear charge of the primray ion before recombination
 * @param E_kin electron-ion collision energy
 * @param nrel principal quantum number up to which levels, rates, and cross sections are treated relativistically
 * @param energies line energies in eV in the resulting line list 
 * @param strengths line stengths in cm2 in the resulting line list 
 * @param widths line widths in eV in the resulting line list
 * @param labels line labels in the resulting line list
 */
static void cascade_spectrum_rel(const RADRATE &hydro, int nmax, double z,
                                 double E_kin, int nrel,
                                 vector<double> &energies,
                                 vector<double> &strengths,
                                 vector<double> &widths,
                                 vector<string> &labels) {
  using hydroconst::hbar_eV_s;

  int ns = nstates_total(nmax, nrel);
  vector<StateInfo> states;
  build_state_info(states, nmax, nrel);

  matrix<double> M(ns, ns, 0.0);

  // Build rate matrix
  for (int i = 0; i < ns; i++) {
    const StateInfo &si = states[i];

    if (si.is_rel) {
      // Relativistic state: use Dirac rates
      double tau = dirac_lifetime(z, si.n, si.kappa, nmax);
      if (tau <= 0.0) continue;
      M(i, i) = -1.0 / tau;

      // Decays to lower relativistic states (E1, E2 and M1 as individual
      // matrix entries)
      for (int j = 0; j < i; j++) {
        const StateInfo &sj = states[j];
        if (!sj.is_rel) continue;
        double tr = dirac_total_transrate(z, si.n, si.kappa, sj.n, sj.kappa);
        if (tr <= 0.0) continue;
        M(j, i) += tr;
      }
    } else {
      // Nonrelativistic state
      double tau = hydro.life(si.n, si.l);
      if (tau <= 0.0) continue;
      M(i, i) = -1.0 / tau;

      // Decays to lower states
      for (int j = 0; j < i; j++) {
        const StateInfo &sj = states[j];
        if (sj.is_rel) {
          // Nonrelativistic upper -> relativistic lower
          // Split nonrel rate by statistical weight of lower j
          if (abs(si.l - sj.l) != 1) continue;
          double tr_nr = hydro.trans(si.n, si.l, sj.n, sj.l);
          if (tr_nr <= 0.0) continue;
          double w = (2.0 * sj.j + 1.0) / (2.0 * (2.0 * sj.l + 1.0));
          M(j, i) += tr_nr * w;
        } else {
          // Nonrelativistic -> nonrelativistic
          if (abs(si.l - sj.l) != 1) continue;
          double br = hydro.branch(si.n, si.l, sj.n, sj.l);
          if (br > 0.0) M(j, i) += br / tau;
        }
      }
    }
  }

  // RR cross sections
  double E_use = (E_kin > 1e-10) ? E_kin : 1e-12;
  vector<double> sigmaRR(ns, 0.0);
  for (int i = 0; i < ns; i++) {
    const StateInfo &s = states[i];
    if (s.is_rel)
      sigmaRR[i] = (2.0 * s.j + 1.0) *
                   dirac_rr_nonneg(dirac_rr_xsec_retarded(E_use, z, s.n, s.kappa),
                                   E_use, s.n, s.kappa);
    else
      sigmaRR[i] = sigmarrqm(E_use, z, s.n, s.l) / E_use;
  }

  // Solve cascade equations (excl. ground state i=0)
  int nexc = ns - 1;
  matrix<double> Mexc(nexc, nexc, 0.0);
  vector<double> rhs(nexc, 0.0);
  for (int i = 0; i < nexc; i++) {
    rhs[i] = -sigmaRR[i + 1];
    for (int j = 0; j < nexc; j++)
      Mexc(i, j) = M(i + 1, j + 1);
  }
  vector<double> x = matrixSolve(Mexc, rhs);

  // Build spectrum
  energies.clear();
  strengths.clear();
  widths.clear();
  labels.clear();

  // Direct RR lines
  for (int i = 0; i < ns; i++) {
    const StateInfo &s = states[i];
    double E_RR;
    if (s.is_rel)
      E_RR = E_kin - dirac_binding_energy(z, s.n, s.kappa);
    else
      E_RR = E_kin - diracEnergy(z, s.n, s.l);
    double strength = sigmaRR[i];
    if (strength <= 0.0) continue;

    char buf[64];
    if (s.is_rel) {
      if (s.l > 0)
        snprintf(buf, sizeof(buf), "(E,%d)/(E,%d)->(%d,%d,%d)", s.l - 1,
                 s.l + 1, s.n, s.l, (int)(2 * s.j + 0.5));
      else
        snprintf(buf, sizeof(buf), "(E,1)->(%d,%d,%d)", s.n, s.l,
                 (int)(2 * s.j + 0.5));
    } else if (s.l > 0)
      snprintf(buf, sizeof(buf), "(E,%d)/(E,%d)->(%d,%d)", s.l - 1, s.l + 1,
               s.n, s.l);
    else
      snprintf(buf, sizeof(buf), "(E,1)->(%d,%d)", s.n, s.l);
    energies.push_back(E_RR);
    strengths.push_back(strength);
    widths.push_back(0.0);
    labels.push_back(buf);
  }

  // Cascade lines
  for (int i = 0; i < ns; i++) {
    const StateInfo &si = states[i];
    if (i == 0) continue; // ground state has no cascade
    int idx_upper = i;
    double pop = x[idx_upper - 1];

    double tau1;
    if (si.is_rel)
      tau1 = dirac_lifetime(z, si.n, si.kappa, nmax);
    else
      tau1 = hydro.life(si.n, si.l);

    for (int j = 0; j < i; j++) {
      const StateInfo &sj = states[j];

      double tr = 0.0;
      double E_photon;
      double inv_tau2 = 0.0;

      if (si.is_rel && sj.is_rel) {
        tr = dirac_total_transrate(z, si.n, si.kappa, sj.n, sj.kappa);
        E_photon = dirac_binding_energy(z, si.n, si.kappa) -
                   dirac_binding_energy(z, sj.n, sj.kappa);
        if (sj.n > 1)
          inv_tau2 = 1.0 / dirac_lifetime(z, sj.n, sj.kappa, nmax);
      } else if (!si.is_rel && sj.is_rel) {
        if (abs(si.l - sj.l) != 1) continue;
        double tr_nr = hydro.trans(si.n, si.l, sj.n, sj.l);
        if (tr_nr <= 0.0) continue;
        double w = (2.0 * sj.j + 1.0) / (2.0 * (2.0 * sj.l + 1.0));
        tr = tr_nr * w;
        E_photon = diracEnergy(z, si.n, si.l) -
                   dirac_binding_energy(z, sj.n, sj.kappa);
        if (sj.n > 1)
          inv_tau2 = 1.0 / dirac_lifetime(z, sj.n, sj.kappa, nmax);
      } else if (!si.is_rel && !sj.is_rel) {
        if (abs(si.l - sj.l) != 1) continue;
        tr = hydro.trans(si.n, si.l, sj.n, sj.l);
        E_photon = diracEnergy(z, si.n, si.l) - diracEnergy(z, sj.n, sj.l);
        if (sj.n > 1)
          inv_tau2 = 1.0 / hydro.life(sj.n, sj.l);
      } else {
        continue; // si rel -> sj nonrel cannot happen
      }

      if (tr <= 0.0) continue;
      double strength = tr * pop;
      double width = hbar_eV_s * (1.0 / tau1 + inv_tau2);

      char buf[64];
      if (si.is_rel)
        snprintf(buf, sizeof(buf), "(%d,%d,%d)->", si.n, si.l,
                 (int)(2 * si.j + 0.5));
      else
        snprintf(buf, sizeof(buf), "(%d,%d)->", si.n, si.l);

      size_t pos = strlen(buf);
      if (sj.is_rel)
        snprintf(buf + pos, sizeof(buf) - pos, "(%d,%d,%d)", sj.n, sj.l,
                 (int)(2 * sj.j + 0.5));
      else
        snprintf(buf + pos, sizeof(buf) - pos, "(%d,%d)", sj.n, sj.l);

      energies.push_back(E_photon);
      strengths.push_back(strength);
      widths.push_back(width);
      labels.push_back(buf);
    }
  }
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief entry into the cacade computation with interactive input
 */
void calc_spectrum() {

  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  printf("\n RR photon spectrum including radiative cascades.\n");
  printf(" Output: Lorentzian line list (energy, strength, width)\n");
  printf(" suitable for convolution with sum_of_peaks.\n\n");

  double z;
  printf("\n Give effective nuclear charge ...................: ");
  scanf("%lf", &z);


  int nmax;
  printf("\n Give maximum principal quantum number n .........: ");
  scanf("%d", &nmax);

  if (nmax > nmaxfactorial / 2.0) {
    nmax = int(nmaxfactorial / 2.0);
    printf("\n nmax set to %4d.\n", nmax);
  }

  int nrel;
  printf("\n Give max relativistic n (=0 for fully nonrel.) ..: ");
  scanf("%d", &nrel);

  if (nrel > 0) {
    printf(" \n Using relativistic (Dirac) rates for n <= %d\n", nrel);
    printf(" Using nonrelativistic rates for n > %d\n", nrel);
  } else {
    printf(" Using nonrelativistic rates for all n\n");
  }

  double E_kin=0.0;
  printf("\n Give electron-ion collision energy (eV) .........: ");
  scanf("%lf", &E_kin);

  double E_min=0.0;
  printf("\n Give minumum photon energy of output (eV)........: ");
  scanf("%lf", &E_min);

  string filenameroot, filename;
  printf("\n Give filename for output (*.RRlines) ............: ");
  cin >> filenameroot;

  printf("\n Now calculating hydrogenic transition rates \n");
  RADRATE hydro(nmax, z, 1.0, 2.0, 2.0);

  int nstates = nstates_total(nmax, nrel);
  printf("\n Number of states: %d\n", nstates);

  vector<double> energies, strengths, widths;
  vector<string> labels;

  cascade_spectrum_rel(hydro, nmax, z, E_kin, nrel, energies, strengths,
                       widths, labels);

  int nlines = (int)energies.size();

  vector<int> idx(nlines);
  for (int i = 0; i < nlines; i++)
    idx[i] = i;
  sort(idx.begin(), idx.end(), [&](int a, int b) {
    return energies[a] < energies[b];
  });

  // we may be not be interested in too small photon energies below E_min
  int nstart=0;
  for (int k = 0; k < nlines; k++) {
    int i = idx[k];
    if (energies[i] > E_min) { 
      nstart = k;
      break;
    }
  }
  filename = filenameroot + ".RRlines";
  ofstream fout(filename);
  
  fout << "#####################################################################"
       << endl;
  fout << "### Photon spectrum resulting from the radiative cascade following"
       << endl;
  fout << "### radiative recombination (RR lines + cascade lines)" << endl;
  fout << "###" << endl;
  fout << "###  hydrocal revision     : " << HYDROCAL_REVISION << endl;
  fout << "###               filename : " << filename << endl;
  fout << "###  number of spec. lines : " << nlines-nstart << endl;
  fout << "###      start date & time : " << asctime(timeinfo);
  fout << "###    nuclear charge Zeff : " << z << endl;
  fout << "###   electron energy (eV) : " << E_kin << endl;
  fout << "###                   nmax : " << nmax << endl;
  fout << "###  relativistic for n<=  : " << nrel << endl;
  fout << "### min photon energy (eV) : " << E_min << endl;
  fout << "###          lines skipped : " << nstart << endl;
  fout << "#####################################################################"
       << endl;
  fout << "#### energy (eV) [1]   cross_section (cm^2) [2]   line-width (eV) [3]"
       << endl;
  fout << "###-----------------------------------------------------------------"
       << endl;

  for (int k = nstart; k < nlines; k++) {
    int i = idx[k];
    fout << "### " << labels[i] << "\n";
    fout << setw(16) << setprecision(8) << scientific << energies[i] << " "
         << setw(16) << setprecision(8) << strengths[i] << " " << setw(16)
         << setprecision(8) << widths[i] << "\n";
  }
  fout.close();
  cout << endl;
  if (E_min > 0.0) {
    cout << " Minimum photon energy for output: " << E_min << " eV,";
    cout << " skipped " << nstart << " lines.";
  }
  cout << endl;
  printf(" Wrote %d lines to %s\n", nlines-nstart, filename.c_str());
  printf(" Use sum_of_peaks (mode 2 = Lorentz) to convolve.\n");
}
