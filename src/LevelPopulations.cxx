/**
 * @file LevelPopulations.cxx
 *
 * @brief Level population dynamics of stored ions
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
 @endverbatim
 *
 */

#include "LevelPopulations.h"
#include "autoso1.h"
#include "buildinfo.h"
#include "hydromath.h"
#include "hydroconst.h"
#include "matrix.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

using namespace std;

namespace {

/** @brief Composes the output file name from an input file name */
string file_prefix(const string &filename) {
  size_t slash = filename.find_last_of("/\\");
  string directory = (slash == string::npos) ? "" : filename.substr(0, slash + 1);
  string base = (slash == string::npos) ? filename : filename.substr(slash + 1);
  size_t dot = base.find_last_of('.');
  if (dot != string::npos && dot > 0)
    base = base.substr(0, dot);
  return directory + base;
}
  
/** @brief Spectroscopic letter belonging to an orbital angular momentum */
string l_symbol(int l) {
  static const char *symbols = "SPDFGHIKLMNOQRTUVWXYZ";
  if (l < 0 || l >= (int)strlen(symbols))
    return "?";
  return string(1, symbols[l]);
}

/** @brief Level label including the configuration, e.g. 2s2 2p6 3s2 3p2 3P1 */
string level_label(const AUTOSO1::Level &level,
                   const vector<string> &configurations) {
  ostringstream label;
  if (level.CF > 0 && (size_t)level.CF <= configurations.size() &&
      !configurations[level.CF - 1].empty())
    label << configurations[level.CF - 1] << ' ';
  label << abs(level.S) << l_symbol(level.L) << 0.5 * level.J;
  return label.str();
}

} // namespace

////////////////////////////////////////////////////////////////////////////
/**
 * @brief Sets up and solves the rate equations 
 *
 * The transition rates are read in from autostructure output.
 * The levels are reordered by decreasing energy. In that order the
 * radiative rate matrix is lower triangular, i.e. its eigenvalues are the
 * negative total outflow rates. This makes it possible to write every
 * population as a sum of exponentials (Bateman solution) and to evaluate
 * the solution cheaply at all requested times. If the data contain upward
 * transitions or degenerate outflow rates, the general matrix-exponential
 * solution from hydromath.cxx is used instead.
 */
int compute_level_populations(const string &filename, const LevelPopParams &params,
			      LevelPopResult &result, string outfilename) {
  AUTOSO1 oic;
  if (!oic.read_oic(filename))
    return 0;

  const vector<AUTOSO1::Level> levels = oic.get_levels();
  const vector<AUTOSO1::RadiativeRate> rates = oic.get_radiative_rates();
  const vector<string> configurations = oic.get_configurations();
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

  // levels to be included in the population table, longest lived first
  vector<size_t> plot_levels;
  for (size_t i : result.lifetime_order) {
    if (result.lifetimes[i] > params.lifetime_min)
      plot_levels.push_back(i);
  }

  // write one combined table: the level data (7 columns) followed by the
  // population data (time, one column per plotted level, and the total)
  ostringstream temperature_string;
  temperature_string << setprecision(6) << params.kB_T_eV;
  const string output_filename =
      outfilename + "T" + temperature_string.str() + ".pop";
  ofstream output_file(output_filename);
  if (!output_file.is_open()) {
    cout << "ERROR:: Cannot write " << output_filename << endl;
    return 0;
  }

  const size_t level_columns = 6;
  const size_t population_columns = 2 + plot_levels.size();
  const size_t columns = level_columns + population_columns;

  time_t rawtime;
  time(&rawtime);
  struct tm *timeinfo = localtime(&rawtime);
  auto metadata = [&](const string &label, const string &value) {
    output_file << "###" << setw(27) << label << " : " << value << "\n";
  };
  auto scientific_string = [](double value) {
    ostringstream value_string;
    value_string << scientific << setprecision(6) << value;
    return value_string.str();
  };
  char start_date[32];
  strftime(start_date, sizeof(start_date), "%Y-%m-%d %H:%M:%S", timeinfo);
  output_file << "###" << string(98, '#') << "\n";
  output_file << "### Level population dynamics from radiative rates\n";
  output_file << "###\n";
  metadata("hydrocal revision", HYDROCAL_REVISION);
  metadata("filename", output_filename);
  metadata("input filename", filename);
  metadata("number of data sets", to_string(max(nlevels, ntimes)));
  metadata("start date & time", start_date);
  metadata("number of levels", to_string(nlevels));
  metadata("plotted levels", to_string(plot_levels.size()));
  metadata("electron temperature (eV)", temperature_string.str());
  metadata("first time point (s)", scientific_string(params.t_start));
  metadata("last time point (s)", scientific_string(params.t_end));
  metadata("number of time points", to_string(ntimes));
  metadata("minimum lifetime (s)", scientific_string(params.lifetime_min));
  output_file << "###" << string(98, '#') << "\n";
  output_file << "### K [1]   LV [2]   T [3]   configuration + label [4]"
                 "   (EK-E1)/eV [5]   lifetime/s [6]   time/s [7]";
  for (size_t k = 0; k < plot_levels.size(); k++)
    output_file << "   " << level_label(levels[plot_levels[k]], configurations) << " ["
                << level_columns + 2 + k << "]";
  output_file << "   total [" << columns << "]\n";
  output_file << "###" << string(columns * 6, '-') << "\n";

  for (size_t row = 0; row < max(nlevels, ntimes); row++) {
    if (row < nlevels) {
      size_t i = result.lifetime_order[row];
      const auto &level = levels[i];
      output_file << level.K << ", " << level.LV << ", " << level.T << ", "
                  << level_label(level, configurations) << ", " << scientific
                  << setprecision(6) << level.energy * hydroconst::Ryd_eV
                  << ", ";
      if (isinf(result.lifetimes[i]))
        output_file << "Infinity";
      else
        output_file << result.lifetimes[i];
    } else {
      for (size_t k = 0; k < level_columns; k++)
        output_file << (k > 0 ? ", --" : "--");
    }

    if (row < ntimes) {
      double total = 0.0;
      for (size_t i = 0; i < nlevels; i++)
        total += result.populations[i][row];
      output_file << ", " << scientific << setprecision(8) << result.times[row];
      for (size_t i : plot_levels)
        output_file << ", " << result.populations[i][row];
      output_file << ", " << total;
    } else {
      for (size_t k = 0; k < population_columns; k++)
        output_file << ", --";
    }
    output_file << "\n";
  }
  output_file.close();

  cout << " Solved rate equations for " << nlevels << " levels on " << ntimes
       << " time points" << endl;
  cout << " Level and population data for " << plot_levels.size()
       << " plotted levels written to " << output_filename << endl;

  return 1;
}

////////////////////////////////////////////////////////////////////////////
/**
 * @brief Interactive front end of the calculation
 */
void LevelPopulations(void) {
  cout << "\n Level population dynamics from radiative rates\n";
  string filename;
  cout << "\n Give name of autostructure file (*.oic or *.o1) ...: ";
  cin >> filename;
  string outfilename = file_prefix(filename);
  LevelPopParams params;
  char answer;
  cout << "\n Change default calculation parameters? (y/n) ......: ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    cout << "\n Give electron temperature in eV (0 for statistical): ";
    cin >> params.kB_T_eV;
    cout << "\n Give first time point in s ........................: ";
    cin >> params.t_start;
    cout << "\n Give last time point in s .........................: ";
    cin >> params.t_end;
    cout << "\n Give number of log-spaced time points .............: ";
    int n_times;
    cin >> n_times;
    if (n_times > 0)
      params.n_times = (size_t)n_times;
    cout << "\n Give mimimum life time in s for level to be plotted: ";
    cin >> params.lifetime_min;
    cout << "\n Give output filename <fn> (<fn>T<temp>.pop) .......: ";
    cin >> outfilename;
  }
  if (params.t_end <= params.t_start) {
    cout << "ERROR:: Last time point must be larger than the first one" << endl;
    return;
  }

  LevelPopResult result;
  string output_prefix;
  compute_level_populations(filename, params, result, outfilename);
}
