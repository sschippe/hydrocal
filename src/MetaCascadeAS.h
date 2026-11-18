/**
 * @file MetaCascadeAS.h
 *
 * @brief Population dynamics from autostructure radiative transition rates
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
 @endverbatim
 *
 */

#pragma once

#include <cstddef>
#include <string>
#include <vector>

/** \brief Parameters of a meta-cascade calculation */
struct MetaCascadeParams {
  double kB_T_eV = 0.0; ///< electron temperature (0 => statistical populations)
  double t_start = 1.0e-10; ///< first time point (s)
  double t_end = 1.0e2;     ///< last time point (s)
  std::size_t n_times = 61; ///< number of logarithmically spaced time points
};

/** \brief Results of a meta-cascade calculation */
struct MetaCascadeResult {
  std::vector<double> times;                    ///< time points (s)
  std::vector<double> lifetimes;                ///< radiative lifetimes (s) per level
  std::vector<std::size_t> lifetime_order;      ///< level indices, longest lived first
  std::vector<std::vector<double>> populations; ///< [level][time], input level order
  double max_norm_deviation = 0.0;              ///< max |sum(P) - 1|
};

/** \brief Runs a meta-cascade calculation from an autostructure *.oic file
 *
 * @param filename  name of the autostructure file (including extension)
 * @param params    calculation parameters
 * @param result    on return, time grid, lifetimes, and level populations
 * @param output_prefix on return, output file name without extension
 * @return 1 on success, 0 on error
 */
int compute_meta_cascade(const std::string &filename,
                         const MetaCascadeParams &params,
                         MetaCascadeResult &result,
                         std::string &output_prefix);

/** \brief Interactive entry point (hydrocal ion-beam menu) */
void MetaCascadeAS(void);
