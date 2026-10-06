/**
 * @file readxsec.h
 *
 * @brief Routines for reading various formats of theory data
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
   @endverbatim
 */

#pragma once

#include <string>
#include <vector>

///////////////////////////////////
/**
 * @enum TheoFormats
 *
 * @brief Formats of theoretical cross section files
 */
enum TheoFormats {
  THEO_FORMAT_UNDEFINED = 0,
  THEO_FORMAT_LORENTZ,
  THEO_FORMAT_AUTOSDR,
  THEO_FORMAT_PINDZOLA,
  THEO_FORMAT_BADNELL,
  THEO_FORMAT_COWAN,
  THEO_FORMAT_GRIFFIN_AS,
  THEO_FORMAT_GRIFFIN_N,
  THEO_FORMAT_GRIFFIN_NL,
  THEO_FORMAT_FRITZSCHE,
  THEO_FORMAT_JAC,
  THEO_FORMAT_STEIH,
  THEO_FORMAT_AUTOSRR,
  THEO_FORMAT_AMARO,
  THEO_FORMAT_ADASADF09
};

void info_DRtheory(void);

int open_DRtheory(std::string filename, int &theo_mode, int &number_of_levels);

void read_DRtheory(std::string filename, int theo_mode, int level_number,
                   std::vector<double> &energy, std::vector<double> &strength,
                   std::vector<double> &wLorentz, int &npts, double &emin,
                   double &emax, double eshift = 0.0);

int open_peak_file(std::string &filename, int &theo_mode,
                   int &number_of_levels);

int read_peak_file(std::string filename, int theo_mode, int level_number,
                   std::vector<double> &energy, std::vector<double> &strength,
                   std::vector<double> &qFano, std::vector<double> &wLorentz,
                   std::vector<double> &wGauss, int npts, double &emin,
                   double &emax);

////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief averages Lorentzian cross section over histogram bins
 */
void bin_lorentzian_peaks(int nL, std::vector<double> eL,
                          std::vector<double> sL, std::vector<double> wL,
                          int nbin, double emin, double bin_width,
                          std::vector<double> &binned_energy,
                          std::vector<double> &binned_xsec);
