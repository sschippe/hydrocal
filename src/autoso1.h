/**
 * @file autoso1.h
 *
 * @brief Processing data from autostructure *.o1 output files
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
 @endverbatim
 *
 */
#pragma once

/** @class AUTOSO1
 * @brief Reads, holds and manipulates a nl-block
 *  of data from an autostructure *.o1 file
 *
 */
#include <fstream>
#include <string>
#include <vector>

void TestAutosO1(void);

class AUTOSO1 {

public:
  /** @brief One level from the NLEVEL block of an autostructure file */
  struct Level {
    int K;      ///< configuration index
    int LV;     ///< level number used in the rate tables
    int T;      ///< term index
    int S;      ///< 2S+1
    int L;      ///< orbital angular momentum (0=S, 1=P, ...)
    int J;      ///< 2J
    int CF;     ///< configuration index
    double energy; ///< (EK-E1)/RY
  };

  /** @brief One radiative transition from the RADIATIVE DATA block */
  struct RadiativeRate {
    int CF1, LV1, W1; ///< upper level (configuration, LV, weight)
    int CF2, LV2, W2; ///< lower level (configuration, LV, weight)
    double rate;      ///< AR*SEC (s^-1), may carry a sign in the file
    double energy;    ///< DEL(RYD)
  };

private:
  int nAA, nLV, nAR;

  std::vector<int> AA_CF1, AA_LV1, AA_W1, AA_CF2, AA_LV2;
  std::vector<double> AA_rate, AA_energy;

  std::vector<int> LV_K, LV_LV, LV_T, LV_S, LV_L, LV_J, LV_CF;
  std::vector<double> LV_energy;

  std::vector<int> AR_CF1, AR_LV1, AR_W1, AR_CF2, AR_LV2, AR_W2;
  std::vector<double> AR_rate, AR_energy;

  std::vector<std::string> configurations; ///< configurations indexed by CF-1

  int nRydberg, lRydberg;
  double Zeff;

  std::fstream fo1;

public:
  int open_file(std::string filename); /**< opens the *.o1 file */

  void close_file(void); /**< closes the *.o1 file */

  int read_nl(int nRyd, int lRyd); /**< reads a block of data
                         belonging to a given nl-Rydberg state */

  int summed_rates(int nminHydro, int CF1); /**< manipulates data
         beloging to a given nl-Rydberg state */

  /** @brief Reads the complete level and radiative-rate data of an
   *         autostructure *.oic (or *.o1) file
   *
   * In contrast to #read_nl this reader does not select an nl-Rydberg
   * state. It reads the first NLEVEL block and all radiative transitions
   * that follow it, which is the layout expected by the meta-cascade
   * calculation.
   *
   * @param filename full file name including its extension
   * @return 1 on success, 0 on error
   */
  int read_oic(const std::string &filename);

  /** @brief Returns the levels read by #read_oic */
  std::vector<Level> get_levels(void) const;

  /** @brief Returns the radiative rates read by #read_oic */
  std::vector<RadiativeRate> get_radiative_rates(void) const;

  /** @brief Returns human-readable configurations indexed by CF-1 */
  std::vector<std::string> get_configurations(void) const;
};
