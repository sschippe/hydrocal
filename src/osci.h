/**
 * @bief Hydrogenic dipole oscillator strengths
 *
 * see Kazem Omidvar and Patricia T. Guimares,
 *     The Astrophysical Journal Suplement Series 73 (1990) 555
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
 @endverbatim
 *
*/

#pragma once

#include <vector>

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief bound-bound hydrogenic oscillator strength for transition n1,l1 -> n2,l2
 *
 * @param n1 principal quantum number of upper bound subshell
 * @param l1 orbital angular moment quantum number of upper bound subshell
 * @param n2 principal quantum number of lower bound subshell
 * @param l2 orbital angular moment quantum number of lower bound subshell
 *
 * @return oscillator strength
 */
double fosciBB(int n1, int l1, int n2, int l2);

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief bound continuum hydrogenic oscillator strength per unit Rydberg at zero electron energy
 *
 * @param n principal quantum number of bound subshell
 * @param l orbital angular moment quantum number of bound subshell
 *
 * @return oscillator strength
 */
double fosciBC(int n, int l);

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief bound-continuum hydrogenic oscillator strength per unit Rydberg
 *
 * @param e continuum energy in Ryd
 * @param n principal quantum number of bound subshell
 * @param l orbital angular moment quantum number of bound subshell
 *
 * @return oscillator strength
 */
double fosciBC(double e, int n, int l);

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief bound-continuum radial matrix elements for all l (Gordon recursion).
 * @param k continuum wavenumber
 * @param n2 principal quantum number of bound shell
 * @param Rbcm on exit Rbcm[l+1] holds the l -> l+1 element
 * @param Rbcp on exit Rbcp[l] holds the l -> l-1 element (l >= 1).
 *
 * The output vectors must be pre-sized to n2+1 elements.
 *
 */
void CalcRbc(double k, int n2, std::vector<double> &Rbcp, std::vector<double> &Rbcm);

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief performs tests of bound-bound oscillator strengths
 */
void testOsciBB(void);

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief performs tests of bound-continuum oscillator strengths
 */
void testOsciBC(void);
