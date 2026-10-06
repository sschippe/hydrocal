// SPDX-License-Identifier: MIT

#pragma once

int calc_alpha(int fselect);

////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief merged-beams rate coefficient from semi-classicale hydrogenic RR cross
 * section
 */
double alpharrscl_cooler(double Erel, double ktpar, double ktperp, double z,
                         int nmin, int nmax, int lmin, int nele);
