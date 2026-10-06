/**
 * @file readadas.h
 *
 * @brief Reads and processes ADASadf09 files
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */
#pragma once

#include "hydroconst.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class ADASadf09 {

public:
  ADASadf09(std::string);

  unsigned int get_npts() const { return ui_npts; }
  unsigned int get_nuccharge() const { return ui_nuccharge; }
  unsigned int get_nlvl() const { return ui_nlvl; }
  unsigned int get_nprnti() const { return ui_nprnti; }
  unsigned int get_nprntf() const { return ui_nprntf; }

  std::vector<double> const &get_temperatures() const { return vd_te; }

  std::vector<double> const &get_alfi(int lvl) const { return vd_alfi[lvl]; }

  std::vector<double> const &get_alft() const { return vd_alft; }

  std::vector<double> const &get_alf_sum() const { return vd_alf_sum; }

  void write_ratecoef();

private:
  unsigned int ui_npts;
  unsigned int ui_nuccharge;
  unsigned int ui_nprnti;
  unsigned int ui_nprntf;
  char c3_isoseq[3];
  char c3_coupling[3];

  double d_bwnp;
  std::vector<unsigned int> vui_indpp;
  std::vector<unsigned int> vui_isp;
  std::vector<unsigned int> vui_ilp;
  std::vector<double> vd_xjp;
  std::vector<double> vd_wnpi;
  std::vector<std::string> vs_cfgp;

  double d_bwnr;
  unsigned int ui_nlvl;
  std::vector<unsigned int> vui_indx;
  std::vector<unsigned int> vui_indp;
  std::vector<unsigned int> vui_is;
  std::vector<unsigned int> vui_il;
  std::vector<double> vd_xj;
  std::vector<double> vd_wnrl;
  std::vector<std::string> vs_cfgl;
  std::vector<double> vd_aalp[5];

  std::vector<double> vd_te;
  std::vector<std::vector<double>> vd_alfi;
  std::vector<std::vector<double>> vd_alff;
  std::vector<double> vd_alft;
  std::vector<double> vd_alf_sum;

  std::string s_filename;
};
