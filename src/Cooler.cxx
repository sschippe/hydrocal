/**
 * @file Cooler.cxx
 *
 * @brief Defines the geometries of the various electron coolers around the
world
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */
#include "Cooler.h"
#include "RRDRxsec.h"
#include "hydroconst.h"
#include "hydromath.h"
#include "kinema.h"
#include "readxsec.h"
#include "buildinfo.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

using namespace std;

/////////////////////////////////////////////////////////////
/**
 * @brief constructor, initializes cooler variables
 *
 * @param id storage ring id, for possible choices see @ref StorageRingIDs
 * @param short_init_flag flags abbreviated initialization
 */
COOLER::COOLER(int id, bool short_init_flag) {
  bool ok = true;

  cooling_energy = 0.0;
  ring_circumference = 0.0;
  max_overlap_length = 0.0;
  nominal_overlap_length = 0.0;
  toroid_sampling_length = 0.0;
  drifttube_halflength = 0.0;
  drifttube_radius = 0.0;
  cathode_radius = 0.0;
  expansion_factor = 0.0;
  beam_area = 0.0;
  beam_radius = 0.0;
  cathode_voltage = 0.0;
  electron_current = 0.0;
  offset_angle = 0.0;
  electron_lab_energy = 0.0;
  electron_lab_voltage = 0.0;
  electron_density_times_beta = 0.0;
  ion_charge_mass_ratio = 0;
  solenoid_length = 0.0;
  angle_flag = false;
  drifttube_flag = false;
  drifttube_scan_flag = false;
  cryring_old_flag = false;
  storage_ring_id = id;
  switch (id) {
  case COOLER::CRY_STORAGE_RING:
    strcpy(storage_ring_string, "CRY");
    CRYinit(short_init_flag);
    break;
  case COOLER::CSR_STORAGE_RING:
    strcpy(storage_ring_string, "CSR");
    CSRinit(short_init_flag);
    break;
  case COOLER::ESR_STORAGE_RING:
    strcpy(storage_ring_string, "ESR");
    ESRinit(short_init_flag);
    break;
  case COOLER::CROSSED_BEAMS:
    strcpy(storage_ring_string, "CB");
    CBinit(short_init_flag);
    break;
  default:
    ok = false;
    storage_ring_id = NO_STORAGE_RING;
  }
  if (ok) {
    beam_radius = cathode_radius * sqrt(expansion_factor);  // in cm^2
    beam_area = hydroconst::pi * beam_radius * beam_radius; // in cm^2
    electron_density_times_beta =
        electron_current /
        (hydroconst::e_As * hydroconst::clight_cm_s * beam_area); // in cm^-3
    cooling_energy = EscFromElab(cathode_voltage, electron_current,
                                 drifttube_radius / beam_radius);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *  @brief passes cooling energy to electron cooler
 *
 * @param Ecool cooling energy in eV
 */
void COOLER::SetCoolingEnergy(double Ecool) {
  cooling_energy = Ecool;
  cooling_voltage = ElabFromEsc(Ecool, electron_current, TubeBeamRatio());
  if (storage_ring_id == CRY_STORAGE_RING) {
    cathode_voltage = cooling_voltage;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief passes the laboratory energy to the electron cooler
 *
 * @param Elabsc laboratory energy in eV (including space charge)
 */
void COOLER::SetLaboratoryEnergy(double Elabsc) {
  electron_lab_energy = Elabsc; // with space charge
  electron_lab_voltage = ElabFromEsc(electron_lab_energy, electron_current,
                                     TubeBeamRatio()); // without space charge
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief initializes CRYRING cooler geometry and other related quantities
 *
 * @param short_init_flag flags abbreviated initialization
 */
void COOLER::CRYinit(bool short_init_flag) {
  const double drifttube_radius_old =
      5.00; // radius of vacuum tube in cooler in cm (before Oct. 2023)
  const double drifttube_radius_new =
      3.80; // radius of  drift tube in cooler in cm (since Oct. 2023)
  const double drifttube_halflength_old =
      10.5; // in cm half length of driftube (two short tubes left and right,
            // befor OCt. 2023)
  const double drifttube_halflength_new =
      45.0; // in cm half length of drifttube (one long tube, after Oct. 2023)
  const double nominal_overlap_length_cathode_scan =
      90.0; // in cm nominal length of overlap between electron beam and ion
            // beam when scanning with cathode voltage
  const double nominal_overlap_length_drifttube_scan =
      70.0; // in cm nominal length of overlap between electron beam and ion
            // beam when scanning with drifttube voltage

  cathode_radius = 0.2;        // cathode radius in cm
  solenoid_length = 90.0;      // length of straight section in cm
  ring_circumference = 5417.0; // in cm
  nominal_overlap_length = nominal_overlap_length_cathode_scan;
  toroid_sampling_length =
      30.0; // sampling length in the Toroid in z direction in cm
  max_overlap_length = solenoid_length + 2.0 * toroid_sampling_length;

  char answer;
  cout << endl << " Intialization of the CRYRING electron cooler";
  cout << endl << " Give cathode voltage (V) ...............................: ";
  cin >> cathode_voltage;
  cooling_voltage = cathode_voltage;
  cout << endl << " Give beam expansion factor .............................: ";
  cin >> expansion_factor;
  cout << endl << " Give electron current (mA) .............................: ";
  cin >> electron_current;
  electron_current *= 0.001; // conversion mA -> A
  if (short_init_flag)
    return;

  cout << endl << " Scanning with cathode or with (old) drift-tube? (c/o/d) : ";
  cin >> answer;
  drifttube_scan_flag = ((answer == 'd') || (answer == 'D'));
  cryring_old_flag = ((answer == 'o') || (answer == 'O'));
  if (cryring_old_flag) {
    drifttube_radius = drifttube_radius_old;
    drifttube_halflength = drifttube_halflength_old;
    cout << endl
         << " Old drift-tube layout (before October 2023), scanning with "
            "cathode";
  } else {
    drifttube_radius = drifttube_radius_new;
    drifttube_halflength = drifttube_halflength_new;
    cout << endl << " New drift-tube layout (after October 2023),";
    if (drifttube_scan_flag) {
      cout << " scanning with drift-tube";
      nominal_overlap_length = nominal_overlap_length_drifttube_scan;
    } else {
      nominal_overlap_length = nominal_overlap_length_cathode_scan;
      cout << " scanning with cathode";
    }
  }
  cout << endl << " Account for spatially varying drift-tube potential?(y/n): ";
  cin >> answer;
  drifttube_flag = ((answer == 'Y') || (answer == 'y'));
  cout << endl << " Account for angular variation of cooler B field? (y/n) .: ";
  cin >> answer;
  angle_flag = ((answer == 'Y') || (answer == 'y'));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief initializes CSR cooler geometry and other related quantities
 *
 * @param short_init_flag flags abbreviated initialization
 */
void COOLER::CSRinit(bool short_init_flag) {
  // CSR values from Leonard Isberner (Email on 2024-10-24)
  ring_circumference = 3512; // in cm
  solenoid_length = 115.2;   // in cm
  nominal_overlap_length = solenoid_length;
  max_overlap_length = solenoid_length; // was left unset (defaulted to 0)
  drifttube_halflength = 43.5; // in cm
  drifttube_radius = 5.0;      // in cm
  cathode_radius = 0.12;       // cathode radius in cm
  expansion_factor = 30.0;
  drifttube_scan_flag =
      true; // CM energy is varied by scanning drifttube potential

  char answer;

  cout << endl << " Initialization of CSR electron cooler" << endl;
  cout << " Give cathode voltage (V) ...............................: ";
  cin >> cathode_voltage;
  cooling_voltage = cathode_voltage;
  cout << " Give beam expansion factor .............................: ";
  cin >> expansion_factor;
  cout << " Give electron current (mA) .............................: ";
  cin >> electron_current;
  electron_current *= 0.001; // conversion mA -> A
  if (short_init_flag)
    return;

  cout << " Account for spatially varying drift-tube potential (y/n): ";
  cin >> answer;
  drifttube_flag = ((answer == 'Y') || (answer == 'y'));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief initializes ESR cooler geometry and other related quantities
 *
 * @param short_init_flag flags abbreviated initialization
 */
void COOLER::ESRinit(bool short_init_flag) {
  const double sampling_start =
      2.0; // in cm (the smallest z-position that should be sampled)
  ring_circumference = 10834; // in cm
  solenoid_length = 250.0;    // in cm
  nominal_overlap_length = solenoid_length;
  max_overlap_length = solenoid_length - 2.0 * sampling_start;
  drifttube_halflength = 97.0; // in cm
  drifttube_radius = 10.0;     // in cm

  cathode_radius =
      2.54; // cathode radius in cm (value from PhD Thesis of Carsten Brandau)
  expansion_factor = 1.0; // electron beam's magnetic expansion factor
  drifttube_scan_flag =
      true; // CM energy is varied by scanning drifttube potential

  char answer;
  cout << endl << " Intialization of ESR electron cooler" << endl;
  cout << " Give cathode voltage (V) ...............................: ";
  cin >> cathode_voltage;
  cooling_voltage = cathode_voltage;
  cout << " Give electron current (mA) .............................: ";
  cin >> electron_current;
  electron_current *= 0.001; // conversion mA -> A
  if (short_init_flag)
    return;

  cout << " Account for spatially varying drifttube potential? (y/n): ";
  cin >> answer;
  drifttube_flag = ((answer == 'Y') || (answer == 'y'));
  cout << " Account for angular variation of cooler B field?   (y/n): ";
  cin >> answer;
  angle_flag = ((answer == 'Y') || (answer == 'y'));
  if (angle_flag) {
    cout << " Give additional offset angle (mrad) ....................: ";
    cin >> offset_angle;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief initializes crossed beams geometry
 *
 * @param short_init_flag flags abbreviated initialization (not used)
 */
void COOLER::CBinit(bool /*short_init_flag*/) {
  electron_current = 0.0; // no space charge effects
  expansion_factor = 1.0;
  cathode_radius = 1.0;
  drifttube_radius = 1.0;
  cout << endl << " Crossed beams geometry" << endl;
  cout << " Give cooling energy (eV) ...............................: ";
  cin >> cooling_energy;
  cooling_voltage =
      cooling_energy; // there are no space charge effects accounted for
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief prints cooler data to output file
 *
 * @param fout file output stream
 */
void COOLER::Print(fstream &fout) {
  if (!fout.is_open())
    return;
  fout << "###           storage ring : " << storage_ring_string << endl;
  fout << "### ring circumference (m) : " << 0.01 * ring_circumference << endl;
  fout << "###   solonoid length (cm) : " << solenoid_length << endl;
  fout << "###    cathode radius (cm) : " << cathode_radius << endl;
  fout << "###       expansion factor : " << expansion_factor << endl;
  fout << "###       beam radius (cm) : " << beam_radius << endl;
  fout << "###  electron current (mA) : " << 1000.0 * electron_current << endl;
  fout << "###    cathode voltage (V) : " << cathode_voltage << endl;
  if (drifttube_scan_flag) {
    fout << "###    energy scanned with : drift-tube" << endl;
  } else {
    fout << "###    energy scanned with : cathode" << endl;
  }
  fout << "### drift-tube length (cm) : " << 2.0 * drifttube_halflength << endl;
  fout << "### drift-tube radius (cm) : " << drifttube_radius << endl;
  if (storage_ring_id == CRY_STORAGE_RING) {
    fout << "###            Publication : Chin. Phys. C 49, 64001 (2025)" << endl;
    fout << "###                    DOI : https://doi.org/10.1088/1674-1137/adbf81"  << endl;
    if (cryring_old_flag) {
      fout << "###       drifttube layout : old (before Oct.2023)" << endl;
      fout << "### max driftube volt. (V) : "
           << 7919 * electron_current / sqrt(cooling_energy) << endl;
    } else {
      fout << "###       drifttube layout : new (after Oct.2023)" << endl;
      fout << "###            Publication : Eur. Phys. J. Spec. Top. (2026)" << endl;
      fout << "###                    DOI : https://doi.org/10.1140/epjs/s11734-026-02357-0"  << endl;
    }
  }
  if (drifttube_flag) {
    fout << "###             drift tube : yes" << endl;
  } else {
    fout << "###             drift tube : no" << endl;
  }
  if (angle_flag) {
    fout << "### nonzero B-field angles : yes" << endl;
    if (storage_ring_id == COOLER::ESR_STORAGE_RING) {
      fout << "###    offset angle (mrad) : " << offset_angle << endl;
    }
  } else {
    fout << "### nonzero B-field angles : no" << endl;
  }
  fout << "###   nominal overlap (cm) : " << nominal_overlap_length << endl;
  fout << "###       max overlap (cm) : " << max_overlap_length << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief plots position dependence of electron beam angle, drifttube potential,
 * and beam deflection
 */
void COOLER::Plot() {
  // output of angle and drifttube voltage as function of z
  fstream fcooler;
  string filename;
  double Elab = (storage_ring_id == CRY_STORAGE_RING)
                    ? 1000.0 + cooling_energy
                    : 1000.0 + cathode_voltage;

  filename = string(storage_ring_string) + "cooler.dat";
  fcooler.open(filename, fstream::out);
  fcooler << "###########################################################"
          << endl;
  fcooler << "### drift-tube potential and angle in " << storage_ring_string
          << " cooler" << endl;
  if (storage_ring_id == CRY_STORAGE_RING) {
    if (cryring_old_flag) {
      fcooler << "###       drifttube layout : old (before Oct.2023)" << endl;
      fcooler << "### max driftube volt. (V) : "
              << 7919 * electron_current / sqrt(Elab) << endl;
    } else {
      fcooler << "###       drifttube layout : new (after Oct.2023)" << endl;
    }
  }
  if (drifttube_flag) {
    fcooler << "### drift-tube potential : yes" << endl;
  } else {
    fcooler << "### drift-tube potential : no" << endl;
  }
  if (angle_flag) {
    fcooler << "###                angle : yes" << endl;
  } else {
    fcooler << "###                angle : no" << endl;
  }
  fcooler << "###            Elab (eV) : " << Elab << endl;
  fcooler << "###   e-beam radius (cm) : " << beam_radius << endl;
  fcooler << "### solenoid length (cm) : " << solenoid_length << endl;
  fcooler << "###    max. overlap (cm) : " << max_overlap_length << endl;
  fcooler << "#################################################################"
             "###################"
          << endl;
  fcooler << "### z position (mm) [1]   Udrift (V) [2]   angle (mrad) [3]    y "
             "deflection (mm) [4]"
          << endl;
  fcooler << "###--------------------------------------------------------------"
             "-------------------"
          << endl;
  for (double z = 0.0; z <= max_overlap_length; z += 0.1) {
    double Udrift_z = CalcDriftTubeVoltage(z);
    double theta_z = 0.0;
    double y_z = 0.0;
    if (BeamAngle(z, theta_z, y_z)) {
      fcooler << 10.0 * (z - 0.5 * max_overlap_length) << ", " << Udrift_z
              << ", " << 1000.0 * theta_z << ", " << 10.0 * y_z << endl;
    }
  }
  fcooler.close();
}

/////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief returns the postion-dependent electron beam angle relative to the
 * cooler axis
 *
 * @param z position along the cooler axis
 * @param angle on exit, the beam angle in rad
 * @param vertical deflection in cm
 *
 * @retrun true, if position is inside allowed range
 * @return false, if position is outside allowed range
 */
bool COOLER::BeamAngle(double z, double &angle, double &y) {
  bool ok = false;
  y = 0.0;
  angle = 0.0;
  switch (storage_ring_id) {
  case COOLER::CRY_STORAGE_RING:
    ok = CRYangle(z, angle, y);
    break;
  case COOLER::CSR_STORAGE_RING:
    ok = true;
    break;
  case COOLER::ESR_STORAGE_RING:
    ok = ESRangle(z, angle);
    break;
  default:
    ok = false;
  }
  return ok;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculates the position-dependent drift tube voltage in the electron
 * cooler
 *
 * @param z position along the cooler axis
 *
 * @return drifttube voltage in V
 */
double COOLER::CalcDriftTubeVoltage(double z) {
  double Udrift = 0.0;
  switch (storage_ring_id) {
  case COOLER::CRY_STORAGE_RING:
    Udrift = CRYdrifttube(z);
    break;
  case COOLER::CSR_STORAGE_RING:
    Udrift = CSRdrifttube(z);
    break;
  case COOLER::ESR_STORAGE_RING:
    Udrift = ESRdrifttube(z);
    break;
  }
  return Udrift;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculates the position-dependent electron density in the electron
 * cooler
 *
 * @param z position along the cooler axis
 *
 * @return electron density in cm^-3
 */
double COOLER::BeamDensity(double z) {
  if (!drifttube_scan_flag) { // scanning with cathode
    cathode_voltage = electron_lab_voltage;
  }
  double Elabscz =
      electron_lab_energy; // z-dependent lab energy with space charge
  if (drifttube_flag) {
    double drifttube_voltage = CalcDriftTubeVoltage(z);
    Elabscz = EscFromElab(cathode_voltage + drifttube_voltage, electron_current,
                          TubeBeamRatio()); // with space charge
  }
  double gamma = 1.0 + Elabscz / hydroconst::mec2_eV;
  double ele_beta = sqrt(1.0 - 1.0 / (gamma * gamma));
  return ele_beta > 0 ? electron_density_times_beta / ele_beta : 0.0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculates position dependent electron lab velocity change in the
 * electron cooler
 *
 * @param z position along the cooler axis
 * @param ion_gamma relativistic gamma factor of the ion (changes in drifttube)
 * @param ele_beta_x x-component of electron velocity at position z
 * @param ele_beta_y y-component of electron velocity at position z
 * @param ele_beta_z z-component of electron velocity at position z
 *
 * @retrun true, if position is inside allowed range
 * @return false, if position is outside allowed range
 */
bool COOLER::BeamBeta(double z, double &ion_gamma, double &/*ele_beta_x*/,
                      double &ele_beta_y, double &ele_beta_z) {
  if (!drifttube_scan_flag) { // scanning with cathode
    cathode_voltage = electron_lab_voltage;
  }
  double Elabscz =
      electron_lab_energy; // z-dependent lab energy with space charge
  if (drifttube_flag) {
    double Udrift = CalcDriftTubeVoltage(z);
    Elabscz = EscFromElab(cathode_voltage + Udrift, electron_current,
                          TubeBeamRatio()); // with space charge
    ion_gamma -= ion_charge_mass_ratio * Udrift;
  }
  double gamma = 1.0 + Elabscz / hydroconst::mec2_eV;
  double ele_beta = sqrt(1.0 - 1.0 / (gamma * gamma));
  ele_beta_z += ele_beta;
  if (angle_flag) { // account for the angular variation of the magnetic guiding
                    // field
    double theta = 0.0;
    double y_z = 0.0; // not used here
    if (!BeamAngle(z, theta, y_z))
      return false; // return with error flag, if z is out of range
    theta -= offset_angle * 0.001; // conversion mrad -> rad
    double ct = cos(theta);
    double st = sin(theta);
    double old_beta_y = ele_beta_y;
    ele_beta_y = old_beta_y * ct - ele_beta_z * st;
    ele_beta_z = old_beta_y * st + ele_beta_z * ct;
  }
  return true; // normal exit
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculates the position-dependent drifttube voltage in the CRYRING
 * electron cooler
 *
 * setup as of before October 2023
 *
 * @param z position along the cooler axis
 *
 * @return drifttube voltage in V
 */
double COOLER::CRYolddrifttube(double z) {
  // The drift tubes in the CRYRING are also referred to as pickups.
  // The voltage on the pickups should compensate the change in space-charge
  // potential that is caused by their narrower diameter as compared to the
  // vacuum tube of the cooler. The narrower diameter decreases the space charge
  // and accelerates the electrons. The additional potential should thus slow
  // down the electrons.
  cathode_voltage =
      electron_lab_voltage; // the cathode potential is used for ramping
  double Udrift_max =
      7919 * electron_current / sqrt(cooling_voltage); // approximate formula
  double Usp_change = 7919 * electron_current /
                      sqrt(electron_lab_voltage); // approximate formula
  double Udrift1 = Udrift_max;
  double zz = z - toroid_sampling_length -
              drifttube_halflength; // zz=0 in the center of the 1st drifttube
  if (zz <= 0.0) {
    Udrift1 *= 0.5 + 0.5 * tanh(1.318 * (zz + drifttube_halflength) /
                                drifttube_radius);
  } else {
    Udrift1 *= 0.5 - 0.5 * tanh(1.318 * (zz - drifttube_halflength) /
                                drifttube_radius);
  }
  if ((zz > -drifttube_halflength) && (zz < drifttube_halflength))
    Udrift1 -= Usp_change;

  double Udrift2 = Udrift_max;
  zz = z - toroid_sampling_length - solenoid_length +
       drifttube_halflength; // zz=0 in the center of the 2nd drifttube
  if (zz <= 0.0) {
    Udrift2 *= 0.5 + 0.5 * tanh(1.318 * (zz + drifttube_halflength) /
                                drifttube_radius);
  } else {
    Udrift2 *= 0.5 - 0.5 * tanh(1.318 * (zz - drifttube_halflength) /
                                drifttube_radius);
  }
  if ((zz > -drifttube_halflength) && (zz < drifttube_halflength))
    Udrift2 -= Usp_change;

  return Udrift1 + Udrift2;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculates the position-dependent drift tube voltage in the CRYRING
 * electron cooler
 *
 * setup as of after October 2023
 *
 * @param z position along the cooler axis
 *
 * @return drifttube voltage in V
 */
double COOLER::CRYdrifttube(double z) {
  if (cryring_old_flag)
    return CRYolddrifttube(z);
  double Udrift = 0.0;
  if (drifttube_scan_flag) {
    Udrift = electron_lab_voltage - cathode_voltage;
  }
  double zz = z - toroid_sampling_length -
              0.5 * solenoid_length; // zz=0 in the center of the cooler
  if (zz <= 0.0) {
    Udrift *= 0.5 + 0.5 * tanh(1.318 * (zz + drifttube_halflength) /
                               drifttube_radius);
  } else {
    Udrift *= 0.5 - 0.5 * tanh(1.318 * (zz - drifttube_halflength) /
                               drifttube_radius);
  }
  return Udrift;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculates the position-dependent drift tube voltage in the CSR
 * electron cooler
 *
 * Formula obtained from Leonard Isberner via email on 2024-10-25
 *
 * @param z position along the cooler axis
 *
 * @return drifttube voltage in V
 */
double COOLER::CSRdrifttube(double z) {
  /*
    Leo: Im CSR ist das Feld der Driftröhren nicht abgeschlossen, sondern geht
    vom (kleineren) Driftrohrradius zum (größeren) Radius der Vakuumkammer.
    Daher ist die tanh-Funktion nicht komplett zutreffend. Aus einer
    Feldsimulation gibt es stattdessen eine Fitformel, die den Potentialverlauf
    auf der Achse beschreibt. Maybe to do: The change of the drift tube radius
    is so far not accounted for in the calculation of the space charge!
  */
  double Udrift = electron_lab_voltage - cathode_voltage;
  if (drifttube_flag) { // account for position-dependent drift tube potential
    double zz =
        z - 0.5 * max_overlap_length; // zz=0 in the center of the cooler
    double zd = fabs(zz) - drifttube_halflength;
    if (zd < 0) {
      Udrift *= 1.0 / (1.0 + exp(zd * 0.4964));
    } else {
      Udrift *= 1.0 / (1.0 + exp(zd * (0.476 - zd * (0.015 - zd * 0.00044))));
    }
  }
  return Udrift;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculates the position-dependent drift tube voltage in the ESR
 * electron cooler
 *
 * @param z position along the cooler axis
 *
 * @return drifttube voltage in V
 */
double COOLER::ESRdrifttube(double z) {
  double Udrift = electron_lab_voltage - cathode_voltage;
  double zz = z - 0.5 * max_overlap_length; // zz=0 in the center of the cooler
  if (zz <= 0.0) {
    Udrift *= 0.5 + 0.5 * tanh(1.318 * (zz + drifttube_halflength) /
                               drifttube_radius);
  } else {
    Udrift *= 0.5 - 0.5 * tanh(1.318 * (zz - drifttube_halflength) /
                               drifttube_radius);
  }
  return Udrift;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculation of the electron beam angle in the CRYRING electron cooler
 *
 * @param z longitudinal coordinate along the cooler axis (in cm)
 * @param theta on exit angle in rad
 * @param y on exit vertical (with respect to cooler axis) deflection of the
 * electron beam
 *
 * @return true, if ok, i.e., always
 */
bool COOLER::CRYangle(double z, double &theta, double &y) {
  // Claude Krantz' +9,3 mrad correction matching Danared's published data
  // [Phys. Scr.48, 405 (1993)]
  double zmm = 10.0 * fabs(z - 0.5 * max_overlap_length); // conversion cm -> mm
  const double z0 = 0;          // the fit formula is  zmm > z0
  const double A1 = 0.30703794; // in rad
  const double z1 = 609.73215;  // in mm
  const double k1 = 52.11124;   // in mm
  const double A2 = 0.22545754; // in rad
  const double z2 = 745.68523;  // in mm
  const double k2 = 59.52064;   // in mm

  if (zmm > z0) {
    double x1 = (z1 - zmm) / k1;
    double x01 = (z1 - z0) / k1;
    double x2 = (z2 - zmm) / k2;
    double x02 = (z2 - z0) / k2;
    double theta0 = A1 / (1.0 + exp(x01)) + A2 / (1.0 + exp(x02));
    theta = A1 / (1.0 + exp(x1)) + A2 / (1.0 + exp(x2)) - theta0;
    // y = integral_z0^z *theta(z') dz'
    double y1 = (zmm - z0) * (1.0 - theta0) +
                k1 * log((1.0 + exp(x1)) / (1.0 + exp(x01)));
    double y2 = (zmm - z0) * (1.0 - theta0) +
                k2 * log((1.0 + exp(x2)) / (1.0 + exp(x02)));
    y = 0.1 * (A1 * y1 + A2 * y2); // in cm
  } else {
    theta = 0.0;
    y = 0.0;
  }
  return (y < beam_radius);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief measured angle in ESR cooler due to imperfections of magnetic guiding
 * field
 *
 * @param zpos position in cooler along z-axis in cm (0 < zpos < 250)
 * @param angle on exit: angle in rad at position zpos
 *
 * @return false, if zpos out of range
 * @return true, if zpos is ok
 */
bool COOLER::ESRangle(double zpos, double &angle) {
  // the angles have been corrected for an overall offset of  0.48 mrad
  const double theta[2501] = {
      0.06216,     0.06145,     0.06075,     0.06004,     0.05934,
      0.05863,     0.05792,     0.05722,     0.05651,     0.05581,
      0.0551,      0.05448,     0.05385,     0.05323,     0.0526,
      0.05198,     0.05135,     0.05072,     0.0501,      0.04947,
      0.04885,     0.04829,     0.04774,     0.04718,     0.04663,
      0.04607,     0.04552,     0.04496,     0.04441,     0.04385,
      0.0433,      0.04281,     0.04232,     0.04183,     0.04134,
      0.04084,     0.04035,     0.03986,     0.03937,     0.03888,
      0.03839,     0.03795,     0.03752,     0.03708,     0.03665,
      0.03621,     0.03577,     0.03534,     0.0349,      0.03447,
      0.03403,     0.03364,     0.03326,     0.03287,     0.03249,
      0.0321,      0.03171,     0.03133,     0.03094,     0.03056,
      0.03017,     0.02983,     0.02948,     0.02914,     0.0288,
      0.02846,     0.02811,     0.02777,     0.02743,     0.02708,
      0.02674,     0.02644,     0.02613,     0.02583,     0.02553,
      0.02522,     0.02492,     0.02462,     0.02432,     0.02401,
      0.02371,     0.02344,     0.02317,     0.0229,      0.02263,
      0.02236,     0.0221,      0.02183,     0.02156,     0.02129,
      0.02102,     0.02078,     0.02054,     0.0203,      0.02006,
      0.01983,     0.01959,     0.01935,     0.01911,     0.01887,
      0.01863,     0.01842,     0.01821,     0.018,       0.01779,
      0.01758,     0.01736,     0.01715,     0.01694,     0.01673,
      0.01652,     0.01633,     0.01614,     0.01596,     0.01577,
      0.01558,     0.01539,     0.0152,      0.01502,     0.01483,
      0.01464,     0.01447,     0.01431,     0.01414,     0.01398,
      0.01381,     0.01364,     0.01348,     0.01331,     0.01315,
      0.01298,     0.01283,     0.01269,     0.01254,     0.01239,
      0.01224,     0.0121,      0.01195,     0.0118,      0.01166,
      0.01151,     0.01138,     0.01125,     0.01112,     0.01099,
      0.01086,     0.01072,     0.01059,     0.01046,     0.01033,
      0.0102,      0.01008,     0.00997,     0.00985,     0.00974,
      0.00962,     0.0095,      0.00939,     0.00927,     0.00916,
      0.00904,     0.00894,     0.00884,     0.00873,     0.00863,
      0.00853,     0.00843,     0.00833,     0.00822,     0.00812,
      0.00802,     0.00793,     0.00784,     0.00775,     0.00766,
      0.00756,     0.00747,     0.00738,     0.00729,     0.0072,
      0.00711,     0.00703,     0.00695,     0.00687,     0.00679,
      0.00671,     0.00662,     0.00654,     0.00646,     0.00638,
      0.0063,      0.00623,     0.00616,     0.00608,     0.00601,
      0.00594,     0.00587,     0.0058,      0.00572,     0.00565,
      0.00558,     0.00552,     0.00545,     0.00539,     0.00533,
      0.00527,     0.0052,      0.00514,     0.00508,     0.00501,
      0.00495,     0.00489,     0.00484,     0.00478,     0.00473,
      0.00467,     0.00461,     0.00456,     0.0045,      0.00445,
      0.00439,     0.00434,     0.00429,     0.00424,     0.00419,
      0.00414,     0.00409,     0.00404,     0.00399,     0.00394,
      0.00389,     0.00385,     0.0038,      0.00376,     0.00371,
      0.00367,     0.00363,     0.00358,     0.00354,     0.00349,
      0.00345,     0.00342,     0.00338,     0.00335,     0.00332,
      0.00328,     0.00325,     0.00322,     0.00319,     0.00315,
      0.00312,     0.00308,     0.00305,     0.00301,     0.00298,
      0.00294,     0.0029,      0.00287,     0.00283,     0.0028,
      0.00276,     0.00273,     0.0027,      0.00266,     0.00263,
      0.0026,      0.00257,     0.00254,     0.0025,      0.00247,
      0.00244,     0.00241,     0.00238,     0.00236,     0.00233,
      0.0023,      0.00227,     0.00224,     0.00222,     0.00219,
      0.00216,     0.00214,     0.00211,     0.00209,     0.00206,
      0.00203,     0.00201,     0.00199,     0.00196,     0.00194,
      0.00191,     0.00189,     0.00187,     0.00184,     0.00182,
      0.0018,      0.00178,     0.00176,     0.00173,     0.00171,
      0.00169,     0.00167,     0.00165,     0.00163,     0.00161,
      0.0016,      0.00158,     0.00156,     0.00154,     0.00152,
      0.0015,      0.00148,     0.00146,     0.00145,     0.00143,
      0.00141,     0.00139,     0.00137,     0.00136,     0.00134,
      0.00132,     0.00131,     0.00129,     0.00128,     0.00126,
      0.00125,     0.00123,     0.00121,     0.0012,      0.00119,
      0.00117,     0.00116,     0.00114,     0.00113,     0.00111,
      0.0011,      0.00109,     0.00107,     0.00106,     0.00104,
      0.00103,     0.00102,     0.00101,     0.000995741, 0.000984321,
      0.000972901, 0.000961481, 0.000950061, 0.000938642, 0.000927222,
      0.000915802, 0.000905462, 0.000895122, 0.000884781, 0.000874441,
      0.000864101, 0.000853761, 0.000843421, 0.00083308,  0.00082274,
      0.0008124,   0.000803041, 0.000793682, 0.000784323, 0.000774964,
      0.000765606, 0.000756247, 0.000746888, 0.000737529, 0.00072817,
      0.000718811, 0.000711031, 0.00070325,  0.00069547,  0.000687689,
      0.000679909, 0.000672129, 0.000664348, 0.000656568, 0.000648787,
      0.000641007, 0.000634271, 0.000627534, 0.000620798, 0.000614062,
      0.000607326, 0.000600589, 0.000593853, 0.000587117, 0.00058038,
      0.000573644, 0.000567285, 0.000560926, 0.000554567, 0.000548208,
      0.000541849, 0.00053549,  0.000529131, 0.000522772, 0.000516413,
      0.000510054, 0.000504092, 0.00049813,  0.000492168, 0.000486206,
      0.000480244, 0.000474282, 0.00046832,  0.000462358, 0.000456396,
      0.000450434, 0.000446974, 0.000443514, 0.000440053, 0.000436593,
      0.000433133, 0.000429673, 0.000426213, 0.000422752, 0.000419292,
      0.000415832, 0.000411735, 0.000407639, 0.000403542, 0.000399445,
      0.000395348, 0.000391252, 0.000387155, 0.000383058, 0.000378962,
      0.000374865, 0.000371337, 0.000367808, 0.00036428,  0.000360751,
      0.000357223, 0.000353695, 0.000350166, 0.000346638, 0.000343109,
      0.000339581, 0.000337572, 0.000335563, 0.000333553, 0.000331544,
      0.000329535, 0.000327526, 0.000325517, 0.000323507, 0.000321498,
      0.000319489, 0.000318455, 0.000317421, 0.000316387, 0.000315353,
      0.000314319, 0.000313284, 0.00031225,  0.000311216, 0.000310182,
      0.000309148, 0.000306429, 0.00030371,  0.000300991, 0.000298272,
      0.000295552, 0.000292833, 0.000290114, 0.000287395, 0.000284676,
      0.000281957, 0.000282048, 0.000282139, 0.00028223,  0.000282321,
      0.000282413, 0.000282504, 0.000282595, 0.000282686, 0.000282777,
      0.000282868, 0.000282918, 0.000282967, 0.000283016, 0.000283066,
      0.000283115, 0.000283165, 0.000283215, 0.000283264, 0.000283314,
      0.000283363, 0.000283326, 0.000283289, 0.000283253, 0.000283216,
      0.000283179, 0.000283142, 0.000283105, 0.000283069, 0.000283032,
      0.000282995, 0.000282928, 0.000282861, 0.000282795, 0.000282728,
      0.000282661, 0.000282594, 0.000282527, 0.000282461, 0.000282394,
      0.000282327, 0.000283347, 0.000284366, 0.000285386, 0.000286405,
      0.000287425, 0.000288444, 0.000289464, 0.000290483, 0.000291503,
      0.000292522, 0.000291682, 0.000290841, 0.000290001, 0.00028916,
      0.00028832,  0.00028748,  0.000286639, 0.000285799, 0.000284958,
      0.000284118, 0.000284701, 0.000285283, 0.000285866, 0.000286449,
      0.000287032, 0.000287614, 0.000288197, 0.00028878,  0.000289362,
      0.000289945, 0.000290553, 0.000291161, 0.000291769, 0.000292377,
      0.000292984, 0.000293592, 0.0002942,   0.000294808, 0.000295416,
      0.000296024, 0.000296657, 0.000297289, 0.000297922, 0.000298555,
      0.000299187, 0.00029982,  0.000300453, 0.000301086, 0.000301718,
      0.000302351, 0.000303532, 0.000304713, 0.000305894, 0.000307075,
      0.000308257, 0.000309438, 0.000310619, 0.0003118,   0.000312981,
      0.000314162, 0.000315679, 0.000317195, 0.000318712, 0.000320228,
      0.000321745, 0.000323261, 0.000324778, 0.000326294, 0.000327811,
      0.000329327, 0.000330194, 0.000331062, 0.000331929, 0.000332797,
      0.000333664, 0.000334531, 0.000335399, 0.000336266, 0.000337134,
      0.000338001, 0.000339075, 0.00034015,  0.000341224, 0.000342299,
      0.000343373, 0.000344447, 0.000345522, 0.000346596, 0.000347671,
      0.000348745, 0.000350017, 0.000351289, 0.000352561, 0.000353833,
      0.000355105, 0.000356377, 0.000357649, 0.000358921, 0.000360193,
      0.000361465, 0.000362366, 0.000363268, 0.000364169, 0.00036507,
      0.000365972, 0.000366873, 0.000367774, 0.000368675, 0.000369577,
      0.000370478, 0.000371106, 0.000371733, 0.000372361, 0.000372989,
      0.000373617, 0.000374244, 0.000374872, 0.0003755,   0.000376127,
      0.000376755, 0.00037699,  0.000377225, 0.00037746,  0.000377695,
      0.000377931, 0.000378166, 0.000378401, 0.000378636, 0.000378871,
      0.000379106, 0.000379303, 0.0003795,   0.000379698, 0.000379895,
      0.000380092, 0.000380289, 0.000380486, 0.000380684, 0.000380881,
      0.000381078, 0.000381259, 0.00038144,  0.00038162,  0.000381801,
      0.000381982, 0.000382163, 0.000382344, 0.000382524, 0.000382705,
      0.000382886, 0.000382704, 0.000382522, 0.00038234,  0.000382158,
      0.000381977, 0.000381795, 0.000381613, 0.000381431, 0.000381249,
      0.000381067, 0.000381154, 0.000381242, 0.000381329, 0.000381416,
      0.000381504, 0.000381591, 0.000381678, 0.000381765, 0.000381853,
      0.00038194,  0.000382784, 0.000383628, 0.000384473, 0.000385317,
      0.000386161, 0.000387005, 0.000387849, 0.000388694, 0.000389538,
      0.000390382, 0.000390399, 0.000390415, 0.000390432, 0.000390448,
      0.000390465, 0.000390481, 0.000390498, 0.000390514, 0.000390531,
      0.000390547, 0.000390991, 0.000391435, 0.000391879, 0.000392323,
      0.000392767, 0.000393212, 0.000393656, 0.0003941,   0.000394544,
      0.000394988, 0.000396004, 0.00039702,  0.000398035, 0.000399051,
      0.000400067, 0.000401083, 0.000402099, 0.000403114, 0.00040413,
      0.000405146, 0.000406156, 0.000407166, 0.000408177, 0.000409187,
      0.000410197, 0.000411207, 0.000412217, 0.000413228, 0.000414238,
      0.000415248, 0.000415928, 0.000416608, 0.000417289, 0.000417969,
      0.000418649, 0.000419329, 0.000420009, 0.00042069,  0.00042137,
      0.00042205,  0.000423291, 0.000424532, 0.000425773, 0.000427014,
      0.000428256, 0.000429497, 0.000430738, 0.000431979, 0.00043322,
      0.000434461, 0.000434986, 0.00043551,  0.000436035, 0.00043656,
      0.000437085, 0.000437609, 0.000438134, 0.000438659, 0.000439183,
      0.000439708, 0.000439665, 0.000439623, 0.00043958,  0.000439538,
      0.000439495, 0.000439452, 0.00043941,  0.000439367, 0.000439325,
      0.000439282, 0.00043967,  0.000440059, 0.000440447, 0.000440836,
      0.000441224, 0.000441612, 0.000442001, 0.000442389, 0.000442778,
      0.000443166, 0.000443157, 0.000443148, 0.000443138, 0.000443129,
      0.00044312,  0.000443111, 0.000443102, 0.000443092, 0.000443083,
      0.000443074, 0.000443173, 0.000443272, 0.000443371, 0.00044347,
      0.000443569, 0.000443668, 0.000443767, 0.000443866, 0.000443965,
      0.000444064, 0.000444131, 0.000444198, 0.000444264, 0.000444331,
      0.000444398, 0.000444465, 0.000444532, 0.000444598, 0.000444665,
      0.000444732, 0.000444783, 0.000444835, 0.000444886, 0.000444938,
      0.000444989, 0.00044504,  0.000445092, 0.000445143, 0.000445195,
      0.000445246, 0.000445545, 0.000445845, 0.000446144, 0.000446443,
      0.000446743, 0.000447042, 0.000447341, 0.00044764,  0.00044794,
      0.000448239, 0.000448746, 0.000449253, 0.000449761, 0.000450268,
      0.000450775, 0.000451282, 0.000451789, 0.000452297, 0.000452804,
      0.000453311, 0.000453799, 0.000454287, 0.000454775, 0.000455263,
      0.000455751, 0.000456239, 0.000456727, 0.000457215, 0.000457703,
      0.000458191, 0.000458927, 0.000459664, 0.000460401, 0.000461137,
      0.000461874, 0.00046261,  0.000463347, 0.000464083, 0.00046482,
      0.000465556, 0.000466107, 0.000466659, 0.00046721,  0.000467762,
      0.000468313, 0.000468864, 0.000469416, 0.000469967, 0.000470519,
      0.00047107,  0.000471506, 0.000471942, 0.000472378, 0.000472814,
      0.00047325,  0.000473685, 0.000474121, 0.000474557, 0.000474993,
      0.000475429, 0.000475697, 0.000475965, 0.000476232, 0.0004765,
      0.000476768, 0.000477036, 0.000477304, 0.000477571, 0.000477839,
      0.000478107, 0.000478097, 0.000478087, 0.000478078, 0.000478068,
      0.000478058, 0.000478048, 0.000478038, 0.000478029, 0.000478019,
      0.000478009, 0.000478112, 0.000478215, 0.000478317, 0.00047842,
      0.000478523, 0.000478626, 0.000478729, 0.000478831, 0.000478934,
      0.000479037, 0.000479305, 0.000479572, 0.00047984,  0.000480108,
      0.000480376, 0.000480643, 0.000480911, 0.000481179, 0.000481446,
      0.000481714, 0.000481755, 0.000481795, 0.000481836, 0.000481876,
      0.000481917, 0.000481957, 0.000481998, 0.000482038, 0.000482079,
      0.000482119, 0.000482206, 0.000482293, 0.000482379, 0.000482466,
      0.000482553, 0.00048264,  0.000482727, 0.000482813, 0.0004829,
      0.000482987, 0.000483052, 0.000483117, 0.000483182, 0.000483247,
      0.000483312, 0.000483376, 0.000483441, 0.000483506, 0.000483571,
      0.000483636, 0.000483845, 0.000484054, 0.000484263, 0.000484472,
      0.000484681, 0.00048489,  0.000485099, 0.000485308, 0.000485517,
      0.000485726, 0.00048624,  0.000486755, 0.000487269, 0.000487784,
      0.000488298, 0.000488812, 0.000489327, 0.000489841, 0.000490356,
      0.00049087,  0.00049129,  0.000491711, 0.000492131, 0.000492552,
      0.000492972, 0.000493392, 0.000493813, 0.000494233, 0.000494654,
      0.000495074, 0.000495396, 0.000495718, 0.000496041, 0.000496363,
      0.000496685, 0.000497007, 0.000497329, 0.000497652, 0.000497974,
      0.000498296, 0.000498734, 0.000499171, 0.000499609, 0.000500046,
      0.000500484, 0.000500921, 0.000501358, 0.000501796, 0.000502234,
      0.000502671, 0.000502891, 0.00050311,  0.00050333,  0.000503549,
      0.000503768, 0.000503988, 0.000504207, 0.000504427, 0.000504646,
      0.000504866, 0.000505106, 0.000505346, 0.000505585, 0.000505825,
      0.000506065, 0.000506305, 0.000506545, 0.000506784, 0.000507024,
      0.000507264, 0.000507319, 0.000507374, 0.000507429, 0.000507484,
      0.000507539, 0.000507594, 0.000507649, 0.000507704, 0.000507759,
      0.000507814, 0.000508189, 0.000508564, 0.00050894,  0.000509315,
      0.00050969,  0.000510065, 0.00051044,  0.000510816, 0.000511191,
      0.000511566, 0.000511841, 0.000512117, 0.000512392, 0.000512667,
      0.000512943, 0.000513218, 0.000513493, 0.000513768, 0.000514044,
      0.000514319, 0.000515365, 0.000516411, 0.000517457, 0.000518503,
      0.000519549, 0.000520595, 0.000521641, 0.000522687, 0.000523733,
      0.000524779, 0.000524736, 0.000524694, 0.000524651, 0.000524609,
      0.000524566, 0.000524523, 0.000524481, 0.000524438, 0.000524396,
      0.000524353, 0.000524421, 0.000524488, 0.000524556, 0.000524623,
      0.000524691, 0.000524759, 0.000524826, 0.000524894, 0.000524961,
      0.000525029, 0.000525076, 0.000525123, 0.000525171, 0.000525218,
      0.000525265, 0.000525312, 0.000525359, 0.000525407, 0.000525454,
      0.000525501, 0.000525978, 0.000526455, 0.000526931, 0.000527408,
      0.000527885, 0.000528362, 0.000528839, 0.000529315, 0.000529792,
      0.000530269, 0.000529511, 0.000528752, 0.000527994, 0.000527236,
      0.000526477, 0.000525719, 0.000524961, 0.000524203, 0.000523444,
      0.000522686, 0.000522719, 0.000522753, 0.000522786, 0.000522819,
      0.000522853, 0.000522886, 0.000522919, 0.000522952, 0.000522986,
      0.000523019, 0.000523272, 0.000523525, 0.000523778, 0.000524031,
      0.000524285, 0.000524538, 0.000524791, 0.000525044, 0.000525297,
      0.00052555,  0.000525733, 0.000525917, 0.0005261,   0.000526283,
      0.000526467, 0.00052665,  0.000526833, 0.000527016, 0.0005272,
      0.000527383, 0.000527496, 0.000527608, 0.000527721, 0.000527833,
      0.000527946, 0.000528059, 0.000528171, 0.000528284, 0.000528396,
      0.000528509, 0.000528299, 0.000528089, 0.000527879, 0.000527669,
      0.000527459, 0.000527249, 0.000527039, 0.000526829, 0.000526619,
      0.000526409, 0.000526598, 0.000526788, 0.000526977, 0.000527166,
      0.000527356, 0.000527545, 0.000527734, 0.000527923, 0.000528113,
      0.000528302, 0.00052832,  0.000528338, 0.000528357, 0.000528375,
      0.000528393, 0.000528411, 0.000528429, 0.000528448, 0.000528466,
      0.000528484, 0.000528522, 0.00052856,  0.000528598, 0.000528636,
      0.000528674, 0.000528713, 0.000528751, 0.000528789, 0.000528827,
      0.000528865, 0.000528603, 0.00052834,  0.000528078, 0.000527815,
      0.000527553, 0.000527291, 0.000527028, 0.000526766, 0.000526503,
      0.000526241, 0.00052686,  0.000527478, 0.000528097, 0.000528715,
      0.000529334, 0.000529953, 0.000530571, 0.00053119,  0.000531808,
      0.000532427, 0.000532679, 0.000532931, 0.000533182, 0.000533434,
      0.000533686, 0.000533938, 0.00053419,  0.000534441, 0.000534693,
      0.000534945, 0.000535001, 0.000535056, 0.000535111, 0.000535167,
      0.000535223, 0.000535278, 0.000535333, 0.000535389, 0.000535445,
      0.0005355,   0.000535288, 0.000535076, 0.000534864, 0.000534652,
      0.00053444,  0.000534228, 0.000534016, 0.000533804, 0.000533592,
      0.00053338,  0.000533182, 0.000532983, 0.000532785, 0.000532586,
      0.000532388, 0.000532189, 0.000531991, 0.000531792, 0.000531594,
      0.000531395, 0.000531162, 0.000530928, 0.000530695, 0.000530462,
      0.000530228, 0.000529995, 0.000529762, 0.000529529, 0.000529295,
      0.000529062, 0.000528601, 0.000528139, 0.000527678, 0.000527217,
      0.000526756, 0.000526294, 0.000525833, 0.000525372, 0.00052491,
      0.000524449, 0.000523872, 0.000523296, 0.000522719, 0.000522142,
      0.000521566, 0.000520989, 0.000520412, 0.000519835, 0.000519259,
      0.000518682, 0.00051864,  0.000518599, 0.000518557, 0.000518516,
      0.000518474, 0.000518432, 0.000518391, 0.000518349, 0.000518308,
      0.000518266, 0.000517943, 0.000517619, 0.000517296, 0.000516973,
      0.00051665,  0.000516326, 0.000516003, 0.00051568,  0.000515356,
      0.000515033, 0.000514652, 0.00051427,  0.000513889, 0.000513507,
      0.000513126, 0.000512744, 0.000512362, 0.000511981, 0.0005116,
      0.000511218, 0.000511405, 0.000511593, 0.00051178,  0.000511967,
      0.000512155, 0.000512342, 0.000512529, 0.000512716, 0.000512904,
      0.000513091, 0.000513156, 0.000513221, 0.000513286, 0.000513351,
      0.000513416, 0.000513481, 0.000513546, 0.000513611, 0.000513676,
      0.000513741, 0.000513685, 0.000513629, 0.000513573, 0.000513517,
      0.000513461, 0.000513405, 0.000513349, 0.000513293, 0.000513237,
      0.000513181, 0.000513419, 0.000513657, 0.000513894, 0.000514132,
      0.00051437,  0.000514608, 0.000514846, 0.000515083, 0.000515321,
      0.000515559, 0.000515797, 0.000516035, 0.000516274, 0.000516512,
      0.00051675,  0.000516988, 0.000517226, 0.000517465, 0.000517703,
      0.000517941, 0.000517974, 0.000518007, 0.00051804,  0.000518073,
      0.000518107, 0.00051814,  0.000518173, 0.000518206, 0.000518239,
      0.000518272, 0.000519111, 0.00051995,  0.000520789, 0.000521628,
      0.000522467, 0.000523305, 0.000524144, 0.000524983, 0.000525822,
      0.000526661, 0.000527174, 0.000527687, 0.0005282,   0.000528713,
      0.000529225, 0.000529738, 0.000530251, 0.000530764, 0.000531277,
      0.00053179,  0.000532255, 0.00053272,  0.000533186, 0.000533651,
      0.000534116, 0.000534581, 0.000535046, 0.000535512, 0.000535977,
      0.000536442, 0.000536616, 0.00053679,  0.000536964, 0.000537138,
      0.000537313, 0.000537487, 0.000537661, 0.000537835, 0.000538009,
      0.000538183, 0.000538034, 0.000537885, 0.000537736, 0.000537587,
      0.000537438, 0.000537288, 0.000537139, 0.00053699,  0.000536841,
      0.000536692, 0.000536108, 0.000535524, 0.00053494,  0.000534356,
      0.000533772, 0.000533187, 0.000532603, 0.000532019, 0.000531435,
      0.000530851, 0.000530599, 0.000530347, 0.000530094, 0.000529842,
      0.00052959,  0.000529338, 0.000529086, 0.000528833, 0.000528581,
      0.000528329, 0.000527972, 0.000527616, 0.00052726,  0.000526903,
      0.000526546, 0.00052619,  0.000525833, 0.000525477, 0.00052512,
      0.000524764, 0.00052466,  0.000524557, 0.000524453, 0.000524349,
      0.000524245, 0.000524142, 0.000524038, 0.000523934, 0.000523831,
      0.000523727, 0.000523255, 0.000522782, 0.00052231,  0.000521838,
      0.000521365, 0.000520893, 0.000520421, 0.000519949, 0.000519476,
      0.000519004, 0.000518694, 0.000518383, 0.000518073, 0.000517762,
      0.000517452, 0.000517142, 0.000516831, 0.000516521, 0.00051621,
      0.0005159,   0.000516428, 0.000516956, 0.000517484, 0.000518012,
      0.00051854,  0.000519067, 0.000519595, 0.000520123, 0.000520651,
      0.000521179, 0.000521662, 0.000522145, 0.000522628, 0.000523111,
      0.000523593, 0.000524076, 0.000524559, 0.000525042, 0.000525525,
      0.000526008, 0.000526765, 0.000527523, 0.00052828,  0.000529038,
      0.000529795, 0.000530552, 0.00053131,  0.000532067, 0.000532825,
      0.000533582, 0.000534644, 0.000535707, 0.000536769, 0.000537831,
      0.000538894, 0.000539956, 0.000541018, 0.00054208,  0.000543143,
      0.000544205, 0.00054402,  0.000543835, 0.00054365,  0.000543465,
      0.00054328,  0.000543095, 0.00054291,  0.000542725, 0.00054254,
      0.000542355, 0.000542067, 0.00054178,  0.000541493, 0.000541205,
      0.000540918, 0.00054063,  0.000540343, 0.000540055, 0.000539768,
      0.00053948,  0.000540022, 0.000540564, 0.000541106, 0.000541648,
      0.000542191, 0.000542733, 0.000543275, 0.000543817, 0.000544359,
      0.000544901, 0.00054451,  0.000544119, 0.000543729, 0.000543338,
      0.000542947, 0.000542556, 0.000542165, 0.000541775, 0.000541384,
      0.000540993, 0.000540454, 0.000539915, 0.000539375, 0.000538836,
      0.000538297, 0.000537758, 0.000537219, 0.000536679, 0.00053614,
      0.000535601, 0.000535888, 0.000536175, 0.000536463, 0.00053675,
      0.000537037, 0.000537324, 0.000537611, 0.000537899, 0.000538186,
      0.000538473, 0.000538137, 0.0005378,   0.000537464, 0.000537128,
      0.000536791, 0.000536455, 0.000536119, 0.000535783, 0.000535446,
      0.00053511,  0.000534523, 0.000533935, 0.000533348, 0.000532761,
      0.000532174, 0.000531586, 0.000530999, 0.000530412, 0.000529824,
      0.000529237, 0.000529089, 0.000528942, 0.000528794, 0.000528646,
      0.000528499, 0.000528351, 0.000528203, 0.000528055, 0.000527908,
      0.00052776,  0.000528237, 0.000528713, 0.00052919,  0.000529666,
      0.000530143, 0.00053062,  0.000531096, 0.000531573, 0.000532049,
      0.000532526, 0.000533038, 0.00053355,  0.000534062, 0.000534574,
      0.000535086, 0.000535598, 0.00053611,  0.000536622, 0.000537134,
      0.000537646, 0.000538005, 0.000538363, 0.000538722, 0.00053908,
      0.000539439, 0.000539798, 0.000540156, 0.000540515, 0.000540873,
      0.000541232, 0.000541044, 0.000540856, 0.000540667, 0.000540479,
      0.000540291, 0.000540103, 0.000539915, 0.000539726, 0.000539538,
      0.00053935,  0.000538807, 0.000538264, 0.00053772,  0.000537177,
      0.000536634, 0.000536091, 0.000535548, 0.000535004, 0.000534461,
      0.000533918, 0.000533905, 0.000533893, 0.00053388,  0.000533867,
      0.000533855, 0.000533842, 0.000533829, 0.000533816, 0.000533804,
      0.000533791, 0.000533486, 0.000533181, 0.000532877, 0.000532572,
      0.000532267, 0.000531962, 0.000531657, 0.000531353, 0.000531048,
      0.000530743, 0.000530166, 0.00052959,  0.000529013, 0.000528437,
      0.00052786,  0.000527283, 0.000526707, 0.00052613,  0.000525554,
      0.000524977, 0.000524357, 0.000523738, 0.000523118, 0.000522498,
      0.000521879, 0.000521259, 0.000520639, 0.000520019, 0.0005194,
      0.00051878,  0.000518603, 0.000518426, 0.000518249, 0.000518072,
      0.000517895, 0.000517718, 0.000517541, 0.000517364, 0.000517187,
      0.00051701,  0.000516063, 0.000515116, 0.000514169, 0.000513222,
      0.000512276, 0.000511329, 0.000510382, 0.000509435, 0.000508488,
      0.000507541, 0.000507323, 0.000507105, 0.000506886, 0.000506668,
      0.00050645,  0.000506232, 0.000506014, 0.000505795, 0.000505577,
      0.000505359, 0.000505332, 0.000505305, 0.000505278, 0.000505251,
      0.000505225, 0.000505198, 0.000505171, 0.000505144, 0.000505117,
      0.00050509,  0.000505341, 0.000505593, 0.000505844, 0.000506096,
      0.000506347, 0.000506598, 0.00050685,  0.000507101, 0.000507353,
      0.000507604, 0.000508071, 0.000508539, 0.000509006, 0.000509473,
      0.000509941, 0.000510408, 0.000510875, 0.000511342, 0.00051181,
      0.000512277, 0.000512978, 0.000513679, 0.000514379, 0.00051508,
      0.000515781, 0.000516482, 0.000517183, 0.000517883, 0.000518584,
      0.000519285, 0.000519377, 0.000519469, 0.000519561, 0.000519653,
      0.000519745, 0.000519837, 0.000519929, 0.000520021, 0.000520113,
      0.000520205, 0.000520173, 0.000520142, 0.00052011,  0.000520078,
      0.000520047, 0.000520015, 0.000519983, 0.000519951, 0.00051992,
      0.000519888, 0.000519905, 0.000519922, 0.000519939, 0.000519956,
      0.000519973, 0.00051999,  0.000520007, 0.000520024, 0.000520041,
      0.000520058, 0.000519915, 0.000519772, 0.000519629, 0.000519486,
      0.000519343, 0.000519201, 0.000519058, 0.000518915, 0.000518772,
      0.000518629, 0.000518054, 0.000517479, 0.000516903, 0.000516328,
      0.000515753, 0.000515178, 0.000514603, 0.000514027, 0.000513452,
      0.000512877, 0.000512773, 0.00051267,  0.000512566, 0.000512462,
      0.000512359, 0.000512255, 0.000512151, 0.000512047, 0.000511944,
      0.00051184,  0.000511824, 0.000511807, 0.000511791, 0.000511775,
      0.000511759, 0.000511742, 0.000511726, 0.00051171,  0.000511693,
      0.000511677, 0.000512367, 0.000513057, 0.000513747, 0.000514437,
      0.000515127, 0.000515818, 0.000516508, 0.000517198, 0.000517888,
      0.000518578, 0.000518656, 0.000518734, 0.000518812, 0.00051889,
      0.000518968, 0.000519046, 0.000519124, 0.000519202, 0.00051928,
      0.000519358, 0.000520128, 0.000520899, 0.000521669, 0.00052244,
      0.00052321,  0.00052398,  0.000524751, 0.000525521, 0.000526292,
      0.000527062, 0.000527607, 0.000528153, 0.000528698, 0.000529244,
      0.000529789, 0.000530334, 0.00053088,  0.000531425, 0.000531971,
      0.000532516, 0.000533411, 0.000534305, 0.0005352,   0.000536094,
      0.000536989, 0.000537883, 0.000538778, 0.000539672, 0.000540566,
      0.000541461, 0.000541601, 0.000541741, 0.000541881, 0.000542021,
      0.000542162, 0.000542302, 0.000542442, 0.000542582, 0.000542722,
      0.000542862, 0.000544118, 0.000545375, 0.000546631, 0.000547887,
      0.000549144, 0.0005504,   0.000551656, 0.000552912, 0.000554169,
      0.000555425, 0.000556221, 0.000557016, 0.000557812, 0.000558608,
      0.000559404, 0.000560199, 0.000560995, 0.000561791, 0.000562586,
      0.000563382, 0.000564182, 0.000564983, 0.000565783, 0.000566583,
      0.000567384, 0.000568184, 0.000568984, 0.000569784, 0.000570585,
      0.000571385, 0.000572218, 0.000573051, 0.000573884, 0.000574717,
      0.00057555,  0.000576384, 0.000577217, 0.00057805,  0.000578883,
      0.000579716, 0.00058035,  0.000580985, 0.000581619, 0.000582254,
      0.000582889, 0.000583523, 0.000584158, 0.000584792, 0.000585427,
      0.000586061, 0.000586728, 0.000587394, 0.000588061, 0.000588728,
      0.000589394, 0.000590061, 0.000590728, 0.000591395, 0.000592061,
      0.000592728, 0.000593476, 0.000594225, 0.000594973, 0.000595722,
      0.00059647,  0.000597218, 0.000597967, 0.000598715, 0.000599464,
      0.000600212, 0.000601466, 0.00060272,  0.000603973, 0.000605227,
      0.000606481, 0.000607735, 0.000608989, 0.000610242, 0.000611496,
      0.00061275,  0.000614047, 0.000615343, 0.00061664,  0.000617936,
      0.000619233, 0.000620529, 0.000621826, 0.000623122, 0.000624419,
      0.000625715, 0.000627924, 0.000630132, 0.000632341, 0.00063455,
      0.000636759, 0.000638967, 0.000641176, 0.000643385, 0.000645593,
      0.000647802, 0.00064976,  0.000651718, 0.000653676, 0.000655634,
      0.000657592, 0.000659549, 0.000661507, 0.000663465, 0.000665423,
      0.000667381, 0.000670074, 0.000672767, 0.00067546,  0.000678153,
      0.000680846, 0.000683539, 0.000686232, 0.000688925, 0.000691618,
      0.000694311, 0.000697682, 0.000701052, 0.000704423, 0.000707794,
      0.000711164, 0.000714535, 0.000717906, 0.000721277, 0.000724647,
      0.000728018, 0.000731152, 0.000734287, 0.000737421, 0.000740556,
      0.000743691, 0.000746825, 0.00074996,  0.000753094, 0.000756229,
      0.000759363, 0.000762827, 0.000766291, 0.000769755, 0.000773219,
      0.000776683, 0.000780147, 0.000783611, 0.000787075, 0.000790539,
      0.000794003, 0.000798409, 0.000802815, 0.000807222, 0.000811628,
      0.000816034, 0.00082044,  0.000824846, 0.000829253, 0.000833659,
      0.000838065, 0.000842554, 0.000847044, 0.000851533, 0.000856022,
      0.000860511, 0.000865001, 0.00086949,  0.000873979, 0.000878469,
      0.000882958, 0.000887472, 0.000891987, 0.000896501, 0.000901016,
      0.00090553,  0.000910044, 0.000914559, 0.000919073, 0.000923588,
      0.000928102, 0.000933822, 0.000939542, 0.000945261, 0.000950981,
      0.000956701, 0.000962421, 0.000968141, 0.00097386,  0.00097958,
      0.0009853,   0.00099177,  0.00099824,  0.001,       0.00101,
      0.00102,     0.00102,     0.00103,     0.00104,     0.00104,
      0.00105,     0.00106,     0.00106,     0.00107,     0.00108,
      0.00109,     0.00109,     0.0011,      0.00111,     0.00111,
      0.00112,     0.00113,     0.00114,     0.00114,     0.00115,
      0.00116,     0.00117,     0.00118,     0.00118,     0.00119,
      0.0012,      0.00121,     0.00122,     0.00123,     0.00124,
      0.00125,     0.00126,     0.00127,     0.00128,     0.00129,
      0.0013,      0.00131,     0.00132,     0.00133,     0.00134,
      0.00135,     0.00137,     0.00138,     0.00139,     0.0014,
      0.00141,     0.00142,     0.00144,     0.00145,     0.00146,
      0.00148,     0.00149,     0.0015,      0.00151,     0.00153,
      0.00154,     0.00155,     0.00157,     0.00158,     0.0016,
      0.00161,     0.00162,     0.00164,     0.00165,     0.00167,
      0.00168,     0.0017,      0.00171,     0.00173,     0.00174,
      0.00176,     0.00177,     0.00179,     0.0018,      0.00182,
      0.00183,     0.00185,     0.00186,     0.00188,     0.0019,
      0.00192,     0.00193,     0.00195,     0.00197,     0.00198,
      0.002,       0.00202,     0.00204,     0.00206,     0.00208,
      0.0021,      0.00212,     0.00214,     0.00216,     0.00218,
      0.0022,      0.00222,     0.00224,     0.00227,     0.00229,
      0.00231,     0.00233,     0.00235,     0.00238,     0.0024,
      0.00242,     0.00244,     0.00247,     0.00249,     0.00252,
      0.00255,     0.00257,     0.0026,      0.00262,     0.00265,
      0.00267,     0.0027,      0.00273,     0.00276,     0.00279,
      0.00281,     0.00284,     0.00287,     0.0029,      0.00293,
      0.00296,     0.00299,     0.00302,     0.00306,     0.00309,
      0.00312,     0.00315,     0.00318,     0.00322,     0.00325,
      0.00328,     0.00332,     0.00335,     0.00339,     0.00343,
      0.00347,     0.0035,      0.00354,     0.00358,     0.00361,
      0.00365,     0.00369,     0.00373,     0.00377,     0.00381,
      0.00385,     0.00389,     0.00393,     0.00397,     0.00401,
      0.00405,     0.00411,     0.00416,     0.00422,     0.00427,
      0.00433,     0.00439,     0.00444,     0.0045,      0.00455,
      0.00461,     0.00466,     0.00471,     0.00477,     0.00482,
      0.00487,     0.00492,     0.00497,     0.00503,     0.00508,
      0.00513,     0.00519,     0.00524,     0.0053,      0.00536,
      0.00541,     0.00547,     0.00553,     0.00559,     0.00564,
      0.0057,      0.00576,     0.00583,     0.00589,     0.00596,
      0.00602,     0.00608,     0.00615,     0.00621,     0.00628,
      0.00634,     0.00641,     0.00648,     0.00655,     0.00662,
      0.00669,     0.00677,     0.00684,     0.00691,     0.00698,
      0.00705,     0.00713,     0.00721,     0.00729,     0.00737,
      0.00745,     0.00752,     0.0076,      0.00768,     0.00776,
      0.00784,     0.00793,     0.00802,     0.00811,     0.0082,
      0.00829,     0.00837,     0.00846,     0.00855,     0.00864,
      0.00873,     0.00883,     0.00892,     0.00902,     0.00912,
      0.00922,     0.00931,     0.00941,     0.00951,     0.0096,
      0.0097,      0.00981,     0.00992,     0.01003,     0.01014,
      0.01025,     0.01035,     0.01046,     0.01057,     0.01068,
      0.01079,     0.01091,     0.01103,     0.01116,     0.01128,
      0.0114,      0.01152,     0.01164,     0.01177,     0.01189,
      0.01201,     0.01214,     0.01228,     0.01241,     0.01255,
      0.01268,     0.01281,     0.01295,     0.01308,     0.01322,
      0.01335,     0.0135,      0.01365,     0.0138,      0.01395,
      0.0141,      0.01425,     0.0144,      0.01455,     0.0147,
      0.01485,     0.01502,     0.01518,     0.01535,     0.01552,
      0.01569,     0.01585,     0.01602,     0.01619,     0.01635,
      0.01652,     0.01671,     0.01689,     0.01708,     0.01726,
      0.01745,     0.01763,     0.01782,     0.018,       0.01819,
      0.01837,     0.01858,     0.01878,     0.01899,     0.0192,
      0.0194,      0.01961,     0.01982,     0.02003,     0.02023,
      0.02044,     0.02067,     0.0209,      0.02113,     0.02136,
      0.02159,     0.02181,     0.02204,     0.02227,     0.0225,
      0.02273,     0.02298,     0.02324,     0.0235,      0.02375,
      0.024,       0.02426,     0.02451,     0.02477,     0.02502,
      0.02528,     0.02556,     0.02585,     0.02613,     0.02642,
      0.0267,      0.02698,     0.02727,     0.02755,     0.02784,
      0.02812,     0.02844,     0.02875,     0.02907,     0.02938,
      0.0297,      0.03002,     0.03033,     0.03065,     0.03096,
      0.03128,     0.03163,     0.03198,     0.03233,     0.03268,
      0.03304,     0.03339,     0.03374,     0.03409,     0.03444,
      0.03479,     0.03518,     0.03557,     0.03596,     0.03635,
      0.03674,     0.03713,     0.03752,     0.03791,     0.0383,
      0.03869,     0.03912,     0.03956,     0.03999,     0.04043,
      0.04086,     0.04129,     0.04173,     0.04216,     0.0426,
      0.04303,     0.04351,     0.044,       0.04448,     0.04496,
      0.04544,     0.04593,     0.04641,     0.04689,     0.04738,
      0.04786,     0.0484,      0.04894,     0.04947,     0.05001,
      0.05055,     0.05109,     0.05163,     0.05216,     0.0527,
      0.05324,     0.05384,     0.05443,     0.05503,     0.05563,
      0.05622,     0.05682,     0.05742,     0.05802,     0.05861,
      0.05921};

  double zmm = 10.0 * zpos; // position in mm
  int i = floor(zmm);
  if ((i < 0) || (i > 2500))
    return false; // out of range

  angle =
      theta[i] + (theta[i + 1] - theta[i]) * (zmm - (double)i); // interpolation
  return true;                                                  // normal exit
}
