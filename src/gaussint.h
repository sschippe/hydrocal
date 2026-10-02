/**
 * @file gaussint.h
 *
 * @brief Provides abscissae and weights for gaussian quadrature
 *
 * @author Stefan Schippers
 * @verbatim
  $Id: gaussint.h 2039 2026-07-20 07:57:32Z iamp $
// SPDX-License-Identifier: MIT
 @endverbatim
 *
 */

#pragma once

#include <vector>

//----------------------------------------------------------------------
/**
 * @brief abscissae and weights for Gauss-Legendre integration
 */
bool gauss_coef(int n, std::vector<double> &x, std::vector<double> &w);

//----------------------------------------------------------------------
/**
 * @brief abscissae and weights for Gauss-Laguerre integration
 */
bool laguer_coef(int n, std::vector<double> &x, std::vector<double> &w);
