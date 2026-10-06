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
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>

using namespace std;

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
