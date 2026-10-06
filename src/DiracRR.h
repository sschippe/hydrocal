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
double dirac_rr_xsec_retarded(double E_el_eV, double z, int n, int kappa);

///////////////////////////////////////////////////////////////////////////////
/**
 * @brief relativistic hydrogenic RR cross section (dipole approximation)
 */
double dirac_rr_xsec_dipole(double E_el_eV, double z, int n, int kappa);

///////////////////////////////////////////////////////////////////////////////
/**
 * @brief clamp a per-vacancy Dirac RR cross section to its physical range
 *
 * dirac_rr_xsec_retarded and dirac_rr_xsec_dipole never return a negative
 * value (0.0 marks a non-existent state or an evaluation that could not be
 * completed).  The aggregation sites weight them by statistical factors
 * (2j+1), so any failure sentinel returned in future would surface in the
 * printed tables as -2, -6, -10, ... rather than as an error.  Any negative
 * value is therefore treated as unavailable (0.0, the convention already used
 * by sigma_IchiharaEichlerRR) and reported once on stderr.  A NaN is left
 * alone: it cannot be mistaken for data.
 *
 * @param sigma    cross section as returned by dirac_rr_xsec_retarded/dipole
 * @param E_el_eV  electron kinetic energy (eV), for the warning
 * @param n        principal quantum number of the captured state
 * @param kappa    relativistic angular quantum number of the captured state
 * @return         sigma if non-negative (or NaN), 0.0 if negative
 */
inline double dirac_rr_nonneg(double sigma, double E_el_eV, int n, int kappa) {
  if (sigma >= 0.0) return sigma;
  if (!(sigma < 0.0)) return sigma; // NaN: propagate, do not mask
  static bool warned = false;
  if (!warned) {
    warned = true;
    fprintf(stderr,
            "WARNING: negative Dirac RR cross section (%.6g cm^2) at "
            "E=%.6g eV, n=%d, kappa=%d; reporting 0 instead\n",
            sigma, E_el_eV, n, kappa);
  }
  return 0.0;
}
