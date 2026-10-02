/**
 * @file Cooler.h
 *
 * @brief Defines the geometries of the various electron coolers around the
world
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: Cooler.h 2039 2026-07-20 07:57:32Z iamp $
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */

#pragma once

#include <fstream>
#include <iostream>

class COOLER {
public:
  /**
   * @enum StorageRingIDs
   *
   * @brief identifies the storage ring in use
   */
  enum StorageRingIDs {
    NO_STORAGE_RING,
    CRY_STORAGE_RING, // CRYRING
    CSR_STORAGE_RING,
    ESR_STORAGE_RING,
    CROSSED_BEAMS
  };

  static int SelectStorageRing(bool crossed_beams_flag = false) {
    std::cout << std::endl << " " << CRY_STORAGE_RING << ") CRYRING";
    std::cout << std::endl << " " << CSR_STORAGE_RING << ") CSR";
    std::cout << std::endl << " " << ESR_STORAGE_RING << ") ESR";
    if (crossed_beams_flag) {
      std::cout << std::endl << " " << CROSSED_BEAMS << ") crossed beams";
    }
    int id = 0;
    std::cout << std::endl
              << " Make a choice ..........................................: ";
    std::cin >> id;
    return id;
  }

  COOLER(int id, bool short_init_flag = false); // constructur
  ~COOLER() { storage_ring_id = NO_STORAGE_RING; } // destructor

  bool Initialized() { return storage_ring_id != NO_STORAGE_RING; }
  void Print(std::fstream &fout);
  void Plot();

  void SetCoolingEnergy(double value);
  void SetLaboratoryEnergy(double value);
  void SetIonChargeMassRatio(double value) { ion_charge_mass_ratio = value; }
  void SetExpansionFactor(double value) { expansion_factor = value; }
  void SetNominalOverlapLength(double value) { nominal_overlap_length = value; }
  double RingCircumference() { return ring_circumference; }
  double ElectronCurrent() { return electron_current; }
  double ExpansionFactor() { return expansion_factor; }
  double NominalOverlapLength() { return nominal_overlap_length; }
  double MaxOverlapLength() { return max_overlap_length; }
  double CoolingVoltage() { return cooling_voltage; }
  double CoolingEnergy() { return cooling_energy; }
  double CathodeRadius() { return cathode_radius; }
  double TubeBeamRatio() {
    return beam_radius > 1E-3 ? drifttube_radius / beam_radius : 0.0;
  }
  double CalcDriftTubeVoltage(double zpos);
  double DriftTubeVoltage() { return electron_lab_voltage - cathode_voltage; }
  double BeamRadius() { return beam_radius; }
  double BeamDensity(double zpos);
  bool BeamAngle(double zpos, double &angle, double &y);
#pragma omp declare simd
  bool BeamBeta(double zpos, double &ion_gamma, double &ele_beta_x,
                double &ele_beta_y, double &ele_beta_z);

private:
  void CBinit(bool);

  void CRYinit(bool);
  bool CRYangle(double zpos, double &angle, double &y);
  double CRYolddrifttube(double zpos);
  double CRYdrifttube(double zpos);

  void CSRinit(bool);
  double CSRdrifttube(double zpos);

  void ESRinit(bool);
  bool ESRangle(double zpos, double &angle);
  double ESRdrifttube(double zpos);

  int storage_ring_id;
  char storage_ring_string[32];
  double cooling_energy;
  double cooling_voltage;
  double ring_circumference;
  double drifttube_halflength;
  double drifttube_radius;
  double cathode_radius;
  double expansion_factor;
  double beam_area;
  double beam_radius;
  double cathode_voltage;
  double electron_lab_energy;
  double electron_lab_voltage;
  double electron_current;
  double offset_angle;
  double electron_density_times_beta;
  double ion_charge_mass_ratio;

  double solenoid_length;
  double toroid_sampling_length;
  double nominal_overlap_length;
  double max_overlap_length;

  bool angle_flag;
  bool drifttube_flag;
  bool drifttube_scan_flag;
  bool cryring_old_flag;

}; // end class COOLER
