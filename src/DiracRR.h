// SPDX-License-Identifier: MIT
/**
 * @file DiracRR.h
 *
 * @brief Relativistic (Dirac) hydrogenic radiative-recombination cross
 *        sections
 *
 * States are specified by (n, kappa); kappa encodes l and j:
 * kappa = -(j+1/2) for j = l+1/2,  kappa = j+1/2 for j = l-1/2.
 *
 */

#pragma once

#include <cstdio>

/** radiative-recombination cross section (cm^2) for the "exact relativistic"
 * calculation of A. Ichihara and J. Eichler, At. Data Nucl. Data Tables 74,
 * 1 (2000), tabulated for all Z = 1..112 and the nine low-lying captured
 * states 1s1/2, 2s1/2, 2p1/2, 2p3/2, 3s1/2, 3p1/2, 3p3/2, 3d3/2, 3d5/2.
 *
 * Natural cubic spline interpolation is performed in log10(E) against
 * log10(sigma*E): near threshold sigma_RR ~ 1/E, so sigma*E is far smoother
 * than sigma itself.  Tabulated points are in barn and are converted to cm²
 * on return.
 *
 * @param E_el_eV free-electron kinetic energy (eV)
 * @param Z       nuclear charge (1..112)
 * @param n       principal quantum number of the captured state
 * @param l       orbital angular momentum of the captured state
 * @param j       total angular momentum (half-integer, e.g. 1.5)
 * @return        sigma_RR (cm²), or 0 if (Z, n, l, j) is not tabulated or
 *                E_el_eV falls outside the energy range for which the paper
 *                prints digits for that state (the columns stop well below
 *                the table's overall energy range, and clamping to the last
 *                printed point would badly overestimate the rapidly-falling
 *                cross section)
 */
double sigma_IchiharaEichlerRR(double E_el_eV, int Z, int n, int l, double j);

///////////////////////////////////////////////////////////////////////////////
/**
 * @brief relativistic hydrogenic RR cross section (fully retarded)
 */
double dirac_rr_xsec_retarded(double E_el_eV, double z, int n, int kappa, double rel_tol, int dj_max, int &dj_conv);
inline double dirac_rr_xsec_retarded(double E_el_eV, double z, int n, int kappa)
{ // use default values for covergence control
  double rel_tol=1.0E-4;
  int dj_max=99;
  int dj_conv;
  return dirac_rr_xsec_retarded(E_el_eV, z, n, kappa, rel_tol, dj_max, dj_conv);
}

///////////////////////////////////////////////////////////////////////////////
/**
 * @brief relativistic hydrogenic RR cross section (dipole approximation)
 */
double dirac_rr_xsec_dipole(double E_el_eV, double z, int n, int kappa);

