/**
 * @file LevelPopulations.h
 *
 * @brief Lvel population dynamics of stored ions
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

/** \brief Parameters of a level population calculation */
struct LevelPopParams {
  double kB_T_eV = 0.0;       ///< electron temperature (0 => statistical populations)
  double t_start = 1.0e-10;   ///< first time point (s)
  double t_end = 1.0e3;       ///< last time point (s)
  std::size_t n_times = 131;  ///< number of logarithmically spaced time points
  double lifetime_min = 1.0;  ///< minimum lifetime (s) for levels in output
};

/** \brief Results of a level population calculation */
struct LevelPopResult {
  std::vector<double> times;                    ///< time points (s)
  std::vector<double> lifetimes;                ///< radiative lifetimes (s) per level
  std::vector<std::size_t> lifetime_order;      ///< level indices, longest lived first
  std::vector<std::vector<double>> populations; ///< [level][time], input level order
  double max_norm_deviation = 0.0;              ///< max |sum(P) - 1|
};

/** \brief Runs a level population calculation
 *
 * @param filename   name of the autostructure file (including extension)
 * @param params     calculation parameters
 * @param result     on return, time grid, lifetimes, and level populations
 * @param oufilename output file name without extension
 * @return 1 on success, 0 on error
 */
int compute_level_populations(const std::string &filename,
			      const LevelPopParams &params,
			      LevelPopResult &result,
			      std::string outfilename);

/** \brief Interactive entry point (hydrocal ion-beam menu) */
void LevelPopulations(void);
