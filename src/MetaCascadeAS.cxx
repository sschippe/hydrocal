/**
 * @file MetaCascadeAS.cxx
 *
 * @brief C++ port of the Meta-Cascade-AS Mathematica notebook
 *
 * The notebook reads radiative transition rates from an autostructure
 * *.oic file, sets up the rate matrix, solves the coupled rate equations
 * for the level populations, and writes level and population tables.
 * Here the same calculations are performed natively; analytic solutions
 * of the triangular rate matrix replace Mathematica's DSolve.
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
 @endverbatim
 *
 */

#include "MetaCascadeAS.h"
#include "autoso1.h"
#include "hydromath.h"
#include "hydroconst.h"
#include "matrix.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

using namespace std;

namespace {

/** \brief Removes leading and trailing whitespace */
string trim(const string &text) {
  size_t first = text.find_first_not_of(" \t\r\n");
  if (first == string::npos)
    return "";
  size_t last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

/** \brief Composes the output file name from an input file name */
string file_prefix(const string &filename) {
  size_t slash = filename.find_last_of("/\\");
  string directory = (slash == string::npos) ? "" : filename.substr(0, slash + 1);
  string base = (slash == string::npos) ? filename : filename.substr(slash + 1);
  size_t dot = base.find_last_of('.');
  if (dot != string::npos && dot > 0)
    base = base.substr(0, dot);
  return directory + base;
}

/** \brief Spectroscopic letter belonging to an orbital angular momentum */
string l_symbol(int l) {
  static const char *symbols = "SPDFGHIKLMNOQRTUVWXYZ";
  if (l < 0 || l >= (int)strlen(symbols))
    return "?";
  return string(1, symbols[l]);
}

/** \brief Level label in the notation used by the notebook, e.g. 3P1 */
string level_label(const AUTOSO1::Level &level) {
  ostringstream label;
  label << abs(level.S) << l_symbol(level.L) << 0.5 * level.J;
  return label.str();
}

/** \brief Reads a value from stdin, using the default on empty input */
string prompt_string(const string &prompt, const string &default_value) {
  cout << prompt << " [" << default_value << "]: ";
  cout.flush();
  string line;
  if (!getline(cin, line))
    return default_value;
  line = trim(line);
  return line.empty() ? default_value : line;
}

/** \brief Reads a double from stdin, using the default on invalid input */
double prompt_double(const string &prompt, double default_value) {
  string value = prompt_string(prompt, "");
  if (value.empty())
    return default_value;
  char *end = nullptr;
  double result = strtod(value.c_str(), &end);
  if (end == value.c_str() || *end != '\0') {
    cout << "  invalid input, using " << default_value << endl;
    return default_value;
  }
  return result;
}

/** \brief Reads an integer from stdin, using the default on invalid input */
int prompt_int(const string &prompt, int default_value) {
  string value = prompt_string(prompt, "");
  if (value.empty())
    return default_value;
  char *end = nullptr;
  long result = strtol(value.c_str(), &end, 10);
  if (end == value.c_str() || *end != '\0' || result <= 0) {
    cout << "  invalid input, using " << default_value << endl;
    return default_value;
  }
  return (int)result;
}

} // namespace

////////////////////////////////////////////////////////////////////////////
/**
 * @brief Sets up and solves the rate equations from an autostructure file
 *
 * The levels are reordered by decreasing energy. In that order the
 * radiative rate matrix is lower triangular, i.e. its eigenvalues are the
 * negative total outflow rates. This makes it possible to write every
 * population as a sum of exponentials (Bateman solution) and to evaluate
 * the solution cheaply at all requested times. If the data contain upward
 * transitions or degenerate outflow rates, the general matrix-exponential
 * solution from hydromath.cxx is used instead.
 */
int compute_meta_cascade(const string &filename, const MetaCascadeParams &params,
                         MetaCascadeResult &result, string &output_prefix) {
  AUTOSO1 oic;
  if (!oic.read_oic(filename))
    return 0;

  const vector<AUTOSO1::Level> levels = oic.get_levels();
  const vector<AUTOSO1::RadiativeRate> rates = oic.get_radiative_rates();
  const size_t nlevels = levels.size();
  if (nlevels == 0 || rates.empty()) {
    cout << "ERROR:: No levels or radiative rates in " << filename << "!"
         << endl;
    return 0;
  }

  // map the LV numbers used in the rate table onto level indices
  int max_lv = 0;
  for (const auto &level : levels)
    max_lv = max(max_lv, level.LV);
  const size_t not_found = numeric_limits<size_t>::max();
  vector<size_t> lv_index(max_lv + 1, not_found);
  for (size_t i = 0; i < nlevels; i++) {
    if (levels[i].LV >= 0 && levels[i].LV <= max_lv)
      lv_index[levels[i].LV] = i;
  }

  // energy order: index 0 is the highest level
  vector<size_t> energy_order(nlevels);
  iota(energy_order.begin(), energy_order.end(), 0);
  sort(energy_order.begin(), energy_order.end(),
       [&](size_t a, size_t b) { return levels[a].energy > levels[b].energy; });
  vector<size_t> position(nlevels);
  for (size_t i = 0; i < nlevels; i++)
    position[energy_order[i]] = i;

  // build the rate matrix: A[final,initial] += |rate|
  matrix<double> A(nlevels, nlevels, 0.0);
  vector<double> outflow(nlevels, 0.0);
  size_t skipped = 0;
  for (const auto &rate : rates) {
    if (rate.LV1 < 0 || rate.LV1 > max_lv || rate.LV2 < 0 ||
        rate.LV2 > max_lv || lv_index[rate.LV1] == not_found ||
        lv_index[rate.LV2] == not_found) {
      skipped++;
      continue;
    }
    size_t initial = position[lv_index[rate.LV1]];
    size_t final = position[lv_index[rate.LV2]];
    double arad = fabs(rate.rate);
    if (arad <= 0.0)
      continue;
    A(final, initial) += arad;
    outflow[initial] += arad;
  }
  for (size_t i = 0; i < nlevels; i++)
    A(i, i) -= outflow[i];
  if (skipped > 0)
    cout << " Warning: skipped " << skipped << " rates with unknown levels"
         << endl;

  // check whether the analytic (triangular) solution can be used
  bool triangular = true;
  for (size_t i = 0; i < nlevels && triangular; i++) {
    for (size_t j = i + 1; j < nlevels; j++) {
      if (A(i, j) != 0.0) {
        triangular = false;
        break;
      }
    }
  }
  bool distinct_rates = true;
  for (size_t i = 0; i < nlevels && distinct_rates; i++) {
    for (size_t j = i + 1; j < nlevels; j++) {
      double scale = max(1.0, max(fabs(outflow[i]), fabs(outflow[j])));
      if (fabs(outflow[i] - outflow[j]) <= 1.0e-10 * scale) {
        distinct_rates = false;
        break;
      }
    }
  }

  // initial populations: statistical weights or Boltzmann distribution
  vector<double> x0(nlevels, 0.0);
  double normalization = 0.0;
  for (size_t i = 0; i < nlevels; i++) {
    double weight = fabs(levels[i].J + 1.0);
    if (params.kB_T_eV > 0.0) {
      double exponent =
          -levels[i].energy * hydroconst::Ryd_eV / params.kB_T_eV;
      exponent = min(exponent, 700.0); // avoid overflow for negative energies
      weight *= exp(exponent);
    }
    x0[position[i]] = weight;
    normalization += weight;
  }
  if (normalization <= 0.0) {
    cout << "ERROR:: Cannot normalize initial populations!" << endl;
    return 0;
  }
  for (auto &weight : x0)
    weight /= normalization;

  // logarithmic time grid
  const size_t ntimes = max<size_t>(1, params.n_times);
  result.times.resize(ntimes);
  double log_start = log10(max(params.t_start, 1.0e-300));
  double log_end = log10(max(params.t_end, params.t_start));
  for (size_t it = 0; it < ntimes; it++) {
    double fraction = (ntimes == 1) ? 0.0 : (double)it / (ntimes - 1);
    result.times[it] = pow(10.0, log_start + fraction * (log_end - log_start));
  }

  // solve dP/dt = A*P in the energy-ordered basis
  vector<vector<double>> ordered_population(nlevels,
                                            vector<double>(ntimes, 0.0));
  if (triangular && distinct_rates) {
    // Bateman coefficients: P_i(t) = sum_{k<=i} c[i][k] exp(-gamma_k t)
    vector<vector<double>> coefficient(nlevels, vector<double>(nlevels, 0.0));
    for (size_t i = 0; i < nlevels; i++) {
      for (size_t j = 0; j < i; j++) {
        double aij = A(i, j);
        if (aij == 0.0)
          continue;
        for (size_t k = 0; k <= j; k++)
          coefficient[i][k] += aij * coefficient[j][k] /
                               (outflow[i] - outflow[k]);
      }
      double already_fed = 0.0;
      for (size_t k = 0; k < i; k++)
        already_fed += coefficient[i][k];
      coefficient[i][i] = x0[i] - already_fed;
    }
    for (size_t it = 0; it < ntimes; it++) {
      for (size_t i = 0; i < nlevels; i++) {
        double population = 0.0;
        for (size_t k = 0; k <= i; k++)
          population += coefficient[i][k] *
                        exp(-outflow[k] * result.times[it]);
        ordered_population[i][it] = population;
      }
    }
  } else {
    cout << " Using matrix-exponential solver (non-triangular or degenerate "
            "rate matrix)"
         << endl;
    for (size_t it = 0; it < ntimes; it++) {
      matrix<double> propagator = matrixExp(A * result.times[it]);
      for (size_t i = 0; i < nlevels; i++) {
        double population = 0.0;
        for (size_t j = 0; j < nlevels; j++)
          population += propagator(i, j) * x0[j];
        ordered_population[i][it] = population;
      }
    }
  }

  // bring lifetimes and populations back into input-level order
  result.lifetimes.resize(nlevels);
  result.populations.assign(nlevels, vector<double>(ntimes, 0.0));
  for (size_t i = 0; i < nlevels; i++) {
    size_t original = energy_order[i];
    result.lifetimes[original] =
        outflow[i] > 0.0 ? 1.0 / outflow[i] : numeric_limits<double>::infinity();
    result.populations[original] = ordered_population[i];
  }

  result.lifetime_order.resize(nlevels);
  iota(result.lifetime_order.begin(), result.lifetime_order.end(), 0);
  stable_sort(result.lifetime_order.begin(), result.lifetime_order.end(),
              [&](size_t a, size_t b) {
                return result.lifetimes[a] > result.lifetimes[b];
              });

  result.max_norm_deviation = 0.0;
  for (size_t it = 0; it < ntimes; it++) {
    double total = 0.0;
    for (size_t i = 0; i < nlevels; i++)
      total += result.populations[i][it];
    result.max_norm_deviation =
        max(result.max_norm_deviation, fabs(total - 1.0));
  }
  if (result.max_norm_deviation > 1.0e-6)
    cout << " Warning: maximum population normalization error "
         << result.max_norm_deviation << endl;

  // write the level table, longest lived first
  output_prefix = file_prefix(filename);
  const string level_filename = output_prefix + "_lev.dat";
  ofstream level_file(level_filename);
  if (level_file.is_open()) {
    level_file << "# K LV T label 2J (EK-E1)/eV lifetime/s\n";
    for (size_t rank = 0; rank < nlevels; rank++) {
      size_t i = result.lifetime_order[rank];
      const auto &level = levels[i];
      level_file << setw(5) << level.K << " " << setw(5) << level.LV << " "
                 << setw(5) << level.T << " " << setw(8) << level_label(level)
                 << " " << setw(5) << level.J << " " << scientific
                 << setprecision(6) << level.energy * hydroconst::Ryd_eV << " ";
      if (isinf(result.lifetimes[i]))
        level_file << "Infinity";
      else
        level_file << result.lifetimes[i];
      level_file << "\n";
    }
    level_file.close();
  } else {
    cout << "ERROR:: Cannot write " << level_filename << endl;
  }

  // write the population table for the longest-lived levels
  ostringstream temperature_string;
  temperature_string << setprecision(6) << params.kB_T_eV;
  const string population_filename = output_prefix + "_" +
                                     to_string(nlevels) + "_" +
                                     temperature_string.str() + ".dat";
  ofstream population_file(population_filename);
  if (population_file.is_open()) {
    const size_t nplot = min<size_t>(5, nlevels);
    vector<size_t> plot_levels(result.lifetime_order.begin(),
                               result.lifetime_order.begin() + nplot);
    population_file << "# time/s";
    for (size_t i : plot_levels)
      population_file << " " << level_label(levels[i]);
    population_file << " total\n";
    for (size_t it = 0; it < ntimes; it++) {
      double total = 0.0;
      for (size_t i = 0; i < nlevels; i++)
        total += result.populations[i][it];
      population_file << scientific << setprecision(8) << result.times[it];
      for (size_t i : plot_levels)
        population_file << " " << result.populations[i][it];
      population_file << " " << total << "\n";
    }
    population_file.close();
  } else {
    cout << "ERROR:: Cannot write " << population_filename << endl;
  }

  cout << " Solved rate equations for " << nlevels << " levels on " << ntimes
       << " time points" << endl;
  cout << " Level data written to " << level_filename << endl;
  cout << " Populations written to " << population_filename << endl;

  return 1;
}

////////////////////////////////////////////////////////////////////////////
/**
 * @brief Interactive front end of the meta-cascade calculation
 */
void MetaCascadeAS(void) {
  cout << "\n Meta-cascade population dynamics from autostructure radiative "
          "rates\n\n";

  // discard the rest of the line containing the menu choice
  cin.ignore(numeric_limits<streamsize>::max(), '\n');

  string filename =
      prompt_string(" Give name of autostructure file (*.oic or *.o1)", "");
  if (filename.empty()) {
    cout << " No file given." << endl;
    return;
  }

  MetaCascadeParams params;
  params.kB_T_eV = prompt_double(
      " Electron temperature at stripping (eV, 0 -> statistical populations)",
      0.0);
  params.t_start = prompt_double(" First time point (s)", params.t_start);
  params.t_end = prompt_double(" Last time point (s)", params.t_end);
  params.n_times = (size_t)prompt_int(" Number of log-spaced time points",
                                      (int)params.n_times);
  if (params.t_end <= params.t_start) {
    cout << "ERROR:: Last time point must be larger than the first one" << endl;
    return;
  }

  MetaCascadeResult result;
  string output_prefix;
  compute_meta_cascade(filename, params, result, output_prefix);
}
