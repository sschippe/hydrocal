/**
 * @file DiracRate.h
 *
 * @brief Relativistic (Dirac) hydrogenic E1, E2, M1 transition rates,
 *        lifetimes and branching ratios
 *
 * States are specified by (n, kappa); kappa encodes l and j:
 * kappa = -(l+1) for j = l+1/2,  kappa = +l for j = l-1/2.
 *
 * $Id: DiracRate.h 2098 2026-08-06 13:21:05Z iamp $
 // SPDX-License-Identifier: MIT
 */

#pragma once

#include <vector>

/** convert orbital and total angular momentum to the Dirac quantum number
 * kappa
 *
 * @param l orbital angular momentum quantum number
 * @param j total angular momentum (half-integer, e.g. 1.5)
 *
 * @return kappa = -(l+1) for j = l+1/2, kappa = +l for j = l-1/2
 */
int kappa_from_lj(int l, double j);

/** convert the Dirac quantum number kappa to l and j
 *
 * @param kappa relativistic angular quantum number
 * @param l     (output) orbital angular momentum quantum number
 * @param j     (output) total angular momentum (half-integer, e.g. 1.5)
 */
void lj_from_kappa(int kappa, int &l, double &j);

/** inverse mapping of dirac_state_index: state index to (n, kappa)
 *
 * @param idx    sequential state index (0-based)
 * @param n      (output) principal quantum number
 * @param kappa  (output) relativistic angular quantum number
 */
void kappa_from_index(int idx, int &n, int &kappa);

/** check whether a state (n, kappa) exists
 *
 * @param n     principal quantum number
 * @param kappa relativistic angular quantum number
 *
 * @return 1 if l < n, i.e. the state is physical, 0 otherwise
 */
bool valid_kappa(int n, int kappa);

/** Dirac (relativistic) binding energy
 *
 * @param z     nuclear charge
 * @param n     principal quantum number
 * @param kappa relativistic angular quantum number
 *
 * @return binding energy in eV (negative)
 */
double dirac_binding_energy(double z, int n, int kappa);

/** Dirac E1 transition rate
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the initial state
 * @param kappa1 angular quantum number of the initial state
 * @param n2     principal quantum number of the final state
 * @param kappa2 angular quantum number of the final state
 *
 * @return the E1 rate in s^-1, 0 if the transition is disallowed
 */
double dirac_transrate(double z, int n1, int kappa1, int n2, int kappa2);

/** Dirac E2 transition rate
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the initial state
 * @param kappa1 angular quantum number of the initial state
 * @param n2     principal quantum number of the final state
 * @param kappa2 angular quantum number of the final state
 *
 * @return the E2 rate in s^-1, 0 if the transition is disallowed
 */
double dirac_e2_rate(double z, int n1, int kappa1, int n2, int kappa2);

/** Dirac M1 transition rate
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the initial state
 * @param kappa1 angular quantum number of the initial state
 * @param n2     principal quantum number of the final state
 * @param kappa2 angular quantum number of the final state
 *
 * @return the M1 rate in s^-1, 0 if the transition is disallowed
 */
double dirac_m1_rate(double z, int n1, int kappa1, int n2, int kappa2);

/** total transition rate including E1, E2, M1
 *
 * Returns 0 when the transition is forbidden by all multipole selection
 * rules.
 *
 * @param z      nuclear charge
 * @param n1     principal quantum number of the initial state
 * @param kappa1 angular quantum number of the initial state
 * @param n2     principal quantum number of the final state
 * @param kappa2 angular quantum number of the final state
 *
 * @return the total rate in s^-1, 0 if the transition is forbidden
 */
double dirac_total_transrate(double z, int n1, int kappa1, int n2, int kappa2);

/** lifetime including E1, E2, and M1 transitions
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
double dirac_lifetime(double z, int n, int kappa, int nmax);

/** branching ratio of a Dirac E1 transition
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
double dirac_branch(double z, int n1, int kappa1, int n2, int kappa2, int nmax);

/** number of Dirac states up to principal quantum number n
 *
 * @param n principal quantum number
 *
 * @return n*n
 */
int ndirac_states(int n);

/** index of the state (n, kappa) in the counting scheme of all Dirac states
 * up to n
 *
 * @param n     principal quantum number
 * @param kappa relativistic angular quantum number
 *
 * @return sequential index (0-based), inverse of kappa_from_index
 */
int dirac_state_index(int n, int kappa);

/** normalisation constant of the Dirac radial functions
 *
 * @param z     nuclear charge
 * @param n     principal quantum number
 * @param kappa relativistic angular quantum number
 *
 * @return the normalisation constant C
 */
double dirac_normalisation(double z, int n, int kappa);

/** build the polynomial coefficients of F(r) and G(r) for a Dirac state
 *
 * F(r) = sum_{k=0}^{nr}  a[k] * r^k
 * G(r) = sum_{k=0}^{nr}  b[k] * r^k
 * where nr = n - |kappa|.
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
void build_FG_coeffs(int kappa, int nr, double gamma, double eps, double two_lambda, double z,
		     std::vector<double> &a, std::vector<double> &b_coeff);
