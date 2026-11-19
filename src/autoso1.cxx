/**
 * @file autoso1.cxx
 *
 * @brief Processing data from autostructure *.o1 output files
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
 @endverbatim
 *
 */

#include "autoso1.h"
#include "fele.h"
#include "hydroconst.h"
#include "lifetime.h"
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <utility>
#include "stdin_guard.h"

using namespace std;

/** @brief Extracts all purely numeric tokens from a line
 *
 * Autostructure output occasionally contains non-numeric tokens such as
 * 'X' in the middle of a data row. They are ignored here.
 */
static void extract_numbers(const string &line, vector<double> &numbers) {
  numbers.clear();
  istringstream sline(line);
  string token;
  while (sline >> token) {
    char *end = nullptr;
    double value = strtod(token.c_str(), &end);
    if (end != token.c_str() && *end == '\0')
      numbers.push_back(value);
  }
}

/** @brief Spectroscopic letter for an orbital angular momentum */
static string orbital_symbol(int l) {
  static const char *symbols = "spdfghiklmnoqrtuvwxyz";
  if (l < 0 || l >= (int)strlen(symbols))
    return "?";
  return string(1, symbols[l]);
}

/** @brief Subshell number encoded in the single-character CFIC label */
static int configuration_index(char label) {
  if (label >= '1' && label <= '9')
    return label - '0';
  if (label >= 'A' && label <= 'Z')
    return 10 + label - 'A';
  if (label >= 'a' && label <= 'z')
    return 36 + label - 'a';
  return -1;
}

/** @brief Parses a fixed-width integer field of a CFIC row */
static bool fixed_int(const string &field, int &value) {
  size_t start = field.find_first_not_of(' ');
  if (start == string::npos)
    return false;
  for (size_t i = start; i < field.size(); i++)
    if (!isdigit((unsigned char)field[i]))
      return false;
  value = atoi(field.c_str());
  return true;
}

int AUTOSO1::open_file(string filename) {
  filename += ".o1";
  fo1.open(filename, fstream::in);
  if (!fo1.is_open()) {
    cout << "ERROR:: File " << filename << " does not exist!" << endl;
    return 0;
  }
  return 1;
}

void AUTOSO1::close_file(void) { fo1.close(); }

int AUTOSO1::read_nl(int nRyd, int lRyd) {
  AA_CF1.clear();
  AA_LV1.clear();
  AA_W1.clear();
  AA_CF2.clear();
  AA_LV2.clear();
  AA_rate.clear();
  AA_energy.clear();
  LV_K.clear();
  LV_LV.clear();
  LV_T.clear();
  LV_S.clear();
  LV_L.clear();
  LV_J.clear();
  LV_CF.clear();
  LV_energy.clear();
  AR_CF1.clear();
  AR_LV1.clear();
  AR_W1.clear();
  AR_CF2.clear();
  AR_LV2.clear();
  AR_W2.clear();
  AR_rate.clear();
  AR_energy.clear();

  string line;
  char Cdummy;
  int nR, lR, nlevel, NucCharge, Nelectron;

  int n2 = 0, n3 = 0, n4 = 0, where = 0;

  //  read first line
  getline(fo1, line);
  while (!fo1.eof()) {
    istringstream sline(line);
    if (line.compare(0, 5, "  NV=") == 0) { // check for next Rydberg electron
      if (where > 1) { // cout << " scan of file finished" << endl;
        break;
      }
      where = 0;
      sline.seekg(5);
      sline >> nR;
      sline.seekg(15);
      sline >> lR;
      getline(fo1, line);
      if ((nR == nRyd) && (lR == lRyd)) {
        where = 1;
        cout << " Rydberg state nl: " << nRyd << " " << lRyd << endl;
        nRydberg = nR;
        lRydberg = lR; // store Rydberg quantum numbers internally
      }
    } else if (line.compare(0, 26, "        I-S            C-S") ==
               0) { // check for autoionization transition data
      if (where > 0) {
        where = 2;
        sline.seekg(66);
        sline >> NucCharge;
        sline.seekg(73);
        sline >> Nelectron;
        cout << "      nuclear charge: " << NucCharge << endl;
        cout << " number of electrons: " << Nelectron << endl;
        Zeff = NucCharge - Nelectron + 1;
      }
      getline(fo1, line);
    } else if (line.compare(0, 10, "   NLEVEL=") == 0) { // check for level data
      if (where > 0) {
        where = 3;
        sline.seekg(10);
        sline >> nlevel;
        cout << " number of levels: " << nlevel << endl;
      }
      getline(fo1, line);
    } else if (line.compare(0, 26, "        I-S            G-S") ==
               0) { // check for radiative transition data
      if (where > 0) {
        where = 4;
      }
      getline(fo1, line);
    } else if (line.find_first_of("1234567890") <
               line.length()) { // if line containes any numbers
      if (where == 2) {         // read autoionization transition data
        int cf1, lv1, w1, cf2, lv2;
        double rate, energy;
        sline >> cf1 >> lv1 >> w1;
        sline >> cf2 >> lv2 >> Cdummy;
        sline >> rate >> energy;
        AA_CF1.push_back(cf1);
        AA_LV1.push_back(lv1);
        AA_W1.push_back(w1);
        AA_CF2.push_back(cf2);
        AA_LV2.push_back(lv2);
        AA_rate.push_back(rate);
        AA_energy.push_back(energy);
        cout << "autoionization transition data: " << cf1 << " " << lv1 << " "
             << w1;
        cout << " " << cf2 << " " << lv2;
        cout << " " << rate << " " << energy << endl;
        n2++;
      } else if (where == 3) { // read level data
        int k, lv, t, s, l, j, cf;
        double energy;
        sline >> k >> lv >> t >> s;
        sline >> l >> j >> cf;
        sline >> energy;
        LV_K.push_back(k);
        LV_LV.push_back(lv);
        LV_T.push_back(t);
        LV_S.push_back(s);
        LV_L.push_back(l);
        LV_J.push_back(j);
        LV_CF.push_back(cf);
        LV_energy.push_back(energy);
        cout << "level data: " << k << " " << lv << " " << t;
        cout << " " << s << " " << l << " ";
        cout << j << " " << cf << " ";
        cout << energy << endl;
        n3++;
      } else if (where == 4) { // read radiative transition data
        int cf1, lv1, w1, cf2, lv2, w2;
        double rate, energy;
        sline >> cf1 >> lv1 >> w1;
        sline >> cf2 >> lv2 >> w2;
        sline >> rate >> energy;
        AR_CF1.push_back(cf1);
        AR_LV1.push_back(lv1);
        AR_W1.push_back(w1);
        AR_CF2.push_back(cf2);
        AR_LV2.push_back(lv2);
        AR_W2.push_back(w2);
        AR_rate.push_back(rate);
        AR_energy.push_back(energy);
        cout << "radiative ransition data: " << cf1 << " " << lv1 << " " << w1;
        cout << " " << cf2 << " " << lv2 << " ";
        cout << w2;
        cout << " " << rate << " " << energy << endl;
        n4++;
      }
    }
    // read next line
    getline(fo1, line);
  } // end while(!fo1.eof())
  nAA = AA_CF1.size();
  if (nAA > 0)
    nAA--;
  nLV = LV_K.size();
  nAR = AR_CF1.size();
  return 1;
}

int AUTOSO1::read_oic(const string &filename) {
  ifstream fin(filename);
  if (!fin.is_open()) {
    cout << "ERROR:: File " << filename << " does not exist!" << endl;
    return 0;
  }

  AA_CF1.clear();
  AA_LV1.clear();
  AA_W1.clear();
  AA_CF2.clear();
  AA_LV2.clear();
  AA_rate.clear();
  AA_energy.clear();
  LV_K.clear();
  LV_LV.clear();
  LV_T.clear();
  LV_S.clear();
  LV_L.clear();
  LV_J.clear();
  LV_CF.clear();
  LV_energy.clear();
  AR_CF1.clear();
  AR_LV1.clear();
  AR_W1.clear();
  AR_CF2.clear();
  AR_LV2.clear();
  AR_W2.clear();
  AR_rate.clear();
  AR_energy.clear();
  configurations.clear();

  string line;
  vector<double> num;

  // The CFIC block precedes NLEVEL. Autostructure writes each occupied
  // subshell as a two-character occupancy followed by a one-character
  // subshell index; all but the first occupancy are offset by 50.
  fin.clear();
  fin.seekg(0);
  while (getline(fin, line)) {
    if (line.find("CFIC") == string::npos)
      continue;
    vector<pair<int, int> > subshells;
    size_t nlpos = line.find("NL");
    if (nlpos != string::npos) {
      vector<double> nl;
      extract_numbers(line.substr(nlpos + 2), nl);
      for (size_t j = 0; j + 1 < nl.size(); j += 2)
        subshells.push_back(make_pair((int)nl[j], (int)nl[j + 1]));
    }
    while (getline(fin, line)) {
      int cf_number, field;
      if (line.size() < 18 || !fixed_int(line.substr(0, 5), cf_number) ||
          !fixed_int(line.substr(5, 5), field) ||
          !isdigit((unsigned char)line[14]) ||
          !fixed_int(line.substr(15, 2), field))
        break;
      ostringstream config;
      bool first = true;
      for (size_t p = 18; p + 2 <= line.size(); p += 3) {
        int occupancy;
        if (!fixed_int(line.substr(p, 2), occupancy))
          break;
        if (!first)
          occupancy -= 50;
        first = false;
        int index = configuration_index(line[p + 2]);
        if (index <= 0 || (size_t)index > subshells.size())
          continue;
        if (!config.str().empty())
          config << ' ';
        config << subshells[index - 1].first
               << orbital_symbol(subshells[index - 1].second) << occupancy;
      }
      configurations.push_back(config.str());
    }
    break;
  }
  fin.clear();
  fin.seekg(0);

  // find the NLEVEL line and extract the number of levels
  int nlevel = 0;
  while (getline(fin, line)) {
    size_t pos = line.find("NLEVEL=");
    if (pos == string::npos)
      continue;
    istringstream sline(line.substr(pos + 7));
    sline >> nlevel;
    break;
  }
  if (nlevel <= 0) {
    cout << "ERROR:: No NLEVEL entry found in " << filename << "!" << endl;
    return 0;
  }

  auto next_numeric_line = [&](size_t nmin) -> bool {
    while (getline(fin, line)) {
      extract_numbers(line, num);
      if (num.size() >= nmin)
        return true;
    }
    return false;
  };

  // read the level table: K LV T 2S+1 L 2J CF (EK-E1)/RY
  for (int i = 0; i < nlevel; i++) {
    if (!next_numeric_line(8)) {
      cout << "ERROR:: Unexpected end of level data in " << filename << "!"
           << endl;
      return 0;
    }
    LV_K.push_back((int)num[0]);
    LV_LV.push_back((int)num[1]);
    LV_T.push_back((int)num[2]);
    LV_S.push_back((int)num[3]);
    LV_L.push_back((int)num[4]);
    LV_J.push_back((int)num[5]);
    LV_CF.push_back((int)num[6]);
    LV_energy.push_back(num[7]);
  }

  // skip to the radiative data block
  bool found_radiative_data = false;
  while (getline(fin, line)) {
    if (line.find("RADIATIVE DATA") != string::npos) {
      found_radiative_data = true;
      break;
    }
  }
  if (!found_radiative_data) {
    cout << " ERROR:: No RADIATIVE DATA block found in " << filename << "!"
         << endl;
    return 0;
  }

  // read all transition rows: CF LV W CF LV W AR*SEC DEL(RYD) E-G(RYD)
  while (next_numeric_line(9)) {
    AR_CF1.push_back((int)num[0]);
    AR_LV1.push_back((int)num[1]);
    AR_W1.push_back((int)num[2]);
    AR_CF2.push_back((int)num[3]);
    AR_LV2.push_back((int)num[4]);
    AR_W2.push_back((int)num[5]);
    AR_rate.push_back(num[6]);
    AR_energy.push_back(num[7]);
  }

  nAA = 0;
  nLV = (int)LV_K.size();
  nAR = (int)AR_CF1.size();
  cout << " Read " << nLV << " levels and " << nAR << " radiative rates from "
       << filename << endl;

  return 1;
}

vector<AUTOSO1::Level> AUTOSO1::get_levels(void) const {
  vector<Level> levels;
  levels.reserve(LV_K.size());
  for (unsigned i = 0; i < LV_K.size(); i++) {
    Level level;
    level.K = LV_K[i];
    level.LV = LV_LV[i];
    level.T = LV_T[i];
    level.S = LV_S[i];
    level.L = LV_L[i];
    level.J = LV_J[i];
    level.CF = LV_CF[i];
    level.energy = LV_energy[i];
    levels.push_back(level);
  }
  return levels;
}

vector<string> AUTOSO1::get_configurations(void) const {
  return configurations;
}

vector<AUTOSO1::RadiativeRate> AUTOSO1::get_radiative_rates(void) const {
  vector<RadiativeRate> rates;
  rates.reserve(AR_CF1.size());
  for (unsigned i = 0; i < AR_CF1.size(); i++) {
    RadiativeRate rate;
    rate.CF1 = AR_CF1[i];
    rate.LV1 = AR_LV1[i];
    rate.W1 = AR_W1[i];
    rate.CF2 = AR_CF2[i];
    rate.LV2 = AR_LV2[i];
    rate.W2 = AR_W2[i];
    rate.rate = AR_rate[i];
    rate.energy = AR_energy[i];
    rates.push_back(rate);
  }
  return rates;
}

int AUTOSO1::summed_rates(int nminHydro, int CF1) {

  // sum of all hydrogenic radiative transition rate from n,l -> n'l' with n'>=
  // nminHydro
  double tau = hydrolife(Zeff, nRydberg, lRydberg, nminHydro);

  std::vector<double> AAsum1(nLV, 0.0);
  std::vector<double> AAsum2(nLV, 0.0);
  std::vector<double> AAsum3(nLV, 0.0);
  std::vector<double> AAenergy1(nLV, 0.0);
  std::vector<double> AAenergy2(nLV, 0.0);
  std::vector<double> AAenergy3(nLV, 0.0);
  std::vector<double> ARsum(nLV, 0.0);
  std::vector<int> AAcount1(nLV, 0);
  std::vector<int> AAcount2(nLV, 0);
  std::vector<int> AAcount3(nLV, 0);
  std::vector<int> ARcount(nLV, 0);
  int LV, i;
  double weight = 0.0;
  for (i = 0; i < nLV; i++) {
    ARsum[i] = 1 / tau;
  }
  for (i = 0; i < nAA; i++) {
    if (AA_CF1[i] == CF1) {
      LV = AA_LV1[i];
      if (((AA_CF2[i] == -20) || (AA_CF2[i] == -17) || (AA_CF2[i] == -14) ||
           (AA_CF2[i] == -11) || (AA_CF2[i] == -8)) &&
          (AA_energy[i] < 0.1) &&
          (AA_rate[i] > 0)) { // rates for dielectronic capture from 3P0
        cout << LV << " " << AA_energy[i] << endl;
        AAsum1[LV - 1] += fabs(AA_rate[i]);
        AAenergy1[LV - 1] += AA_energy[i];
        AAcount1[LV - 1]++;
      } else if (((AA_CF2[i] == -19) || (AA_CF2[i] == -16) ||
                  (AA_CF2[i] == -13) || (AA_CF2[i] == -10) ||
                  (AA_CF2[i] == -7)) &&
                 (AA_rate[i] > 0)) { // rates for autoionization to the ground
                                     // state cout << AA_energy[i] << endl;
        AAsum2[LV - 1] += AA_rate[i];
        AAenergy2[LV - 1] += AA_energy[i];
        AAcount2[LV - 1]++;
      } else {
        AAsum3[LV - 1] += AA_rate[i];
        AAenergy3[LV - 1] += AA_energy[i];
        AAcount3[LV - 1]++;
      }
    }
  }
  for (i = 0; i < nAR; i++) {
    if ((AR_CF1[i] == CF1) && (AR_rate[i] > 0)) {
      LV = AA_LV1[i];
      ARsum[LV - 1] += AR_rate[i];
      ARcount[LV - 1]++;
    }
  }
  cout << "   # 2J+1       Aa3P0   EDR3P0 ";
  cout << "       Aa1S0     E1S0         Arad";
  cout << "          EDR3P0      SDR      SRE  SRE/SDR  alphaDR" << endl;
  for (i = 0; i < nLV; i++) {
    if ((AAsum1[i] > 0.0)) {
      for (int j = 0; j < nLV; j++) {
        if (LV_LV[j] == i)
          weight = LV_J[j] + 1.0;
      }
      cout.width(4);
      cout.precision(3);
      cout << i << " ";
      cout.width(3);
      cout << weight << " ";
      cout.width(3);
      cout << AAcount1[i] << " ";
      cout.width(8);
      cout << AAsum1[i] << " ";
      cout.width(8);
      cout << hydroconst::Ryd_eV * AAenergy1[i] / AAcount1[i] << " ";
      cout.width(3);
      cout << AAcount2[i] << " ";
      cout.width(8);
      cout << AAsum2[i] << " ";
      cout.width(8);
      cout << hydroconst::Ryd_eV * AAenergy2[i] / AAcount2[i] << " ";
      //            cout.width(3);
      //            cout << AAcount3[i] << " ";
      //            cout.width(8);
      //            cout << AAsum3[i] << " ";
      //            cout.width(8);
      //            cout << hydroconst::Ryd_eV*AAenergy3[i]/AAcount3[i] << " ";
      cout.width(3);
      cout << ARcount[i] << " ";
      cout.width(8);
      cout << ARsum[i] << "        ";
      double vel = sqrt(2 * AAenergy1[i] / hydroconst::mec2_eV) *
                   hydroconst::clight_cm_s;
      ;
      double Edelta =
          4 * sqrt(0.69315 * 5E-5 * AAenergy1[i]); // kTpar = 50 mueV
      double sum = AAsum1[i] + AAsum2[i] + ARsum[i];
      double DRstrength =
          4.95e-30 * weight * AAsum1[i] * ARsum[i] / sum / 2 / AAenergy1[i];
      double REstrength =
          4.95e-30 * weight * AAsum1[i] * AAsum2[i] / sum / 2 / AAenergy1[i];
      double DRalpha =
          DRstrength * vel * 0.939 / Edelta; // maximum of rate coefficient
      cout.width(8);
      cout << hydroconst::Ryd_eV * AAenergy1[i] / AAcount1[i] << " ";
      cout.width(8);
      cout << DRstrength << " ";
      cout.width(8);
      cout << REstrength << " ";
      cout.width(8);
      cout << REstrength / DRstrength << " ";
      cout.width(8);
      cout << DRalpha << endl;
    }
  }

  return 0;
}

void TestAutosO1(void) {
  AUTOSO1 AO1;
  int nRyd, lRyd;

  string filename;

  cout << "Give name of autostructure file (*.o1): ";
  cin >> filename;
  AO1.open_file(filename);

  cout << "Give Rydberg n and l: ";
  cin >> nRyd >> lRyd;
  AO1.read_nl(nRyd, lRyd);
  AO1.summed_rates(3, 5);
  AO1.close_file();
}
