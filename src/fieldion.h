/**
 * @file fieldion.h
 *
 * @brief field ionization survival probabilities
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */

#pragma once

#include "hydroconst.h"
#include <vector>
double const Fau =
    hydroconst::au_F_V_cm;             // atomic unit of field strength in V/cm
double const Tau = hydroconst::au_t_s; // atomic unit of time in s

void testFI(void);

void calc_survival(double dt, double f, double z, int nmax,
                   std::vector<double> &psurv, int &n_one, int &n_zero);
