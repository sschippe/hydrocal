/**
 * @file hydroconst.h
 *
 * @brief mathematical and physical constants
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: hydroconst.h 2039 2026-07-20 07:57:32Z iamp $
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */

#pragma once

namespace hydroconst {

constexpr double pi = 3.14159265358979324;
constexpr double sqrtpi = 1.77245385090551603;
constexpr double sqrtln2 = 0.83255461115769776;
constexpr double clight_nm_fs = 299.792458; ///< speed of light in nm/fs (CODATA 2022)
constexpr double clight_m_s = 299792458.0; ///< speed of light in m/s (CODATA 2022)
constexpr double clight_cm_s = 29979245800.0; ///< speed of light in cm/s (CODATA 2022)
constexpr double alpha = 7.2973525643E-3; ///< fine structure constant (CODATA 2022)
constexpr double e_As = 1.60217634E-19; ///< elementary charge in As (CODATA 2022)
constexpr double Ryd_eV = 13.605693122990; ///< Rydberg constant in eV (CODATA 2022)
constexpr double mec2_eV = 510998.9569; ///< electron rest mass in eV (CODATA 2022)
constexpr double mec2_MeV = 0.5109989569; ///< electron rest mass in MeV (CODATA 2022)
constexpr double muc2_eV = 931494103.72; ///< atomic mass unit in eV (CODATA 2022)
constexpr double muc2_MeV = 931.49410372; ///< atomic mass unit in MeV (CODATA 2022)
constexpr double mpc2_eV = 938272089.43; ///< proton mass in eV (CODATA 2022)
constexpr double mnc2_eV = 939565421.94; ///< neutron mass in eV (CODATA 2022)
constexpr double malphac2_eV = 3727379411.8; ///< alpha particle mass in eV (CODATA 2022)
constexpr double mu_kg = 1.66053906892E-27; ///< atomic mass unit in kg (CODATA 2022)
constexpr double hc_eV_nm = 1239.841984;      ///< h*c in eV nm (CODATA 2022)
constexpr double hbarc_eV_nm = 197.3269804;   ///< hbar*c in eV nm (CODATA 2022)
constexpr double h_eV_s = 4.135667696E-15; ///< h in eV s (CODATA 2022)
constexpr double hbar_eV_s = 6.582119569E-16; ///< hbar in eV s (CODATA 2022)
constexpr double kB_J_K = 1.380649E-23; ///< Boltzmann constant in J/K (CODATA 2022)
constexpr double kB_eV_K = 8.617333262E-5; ///< Boltzmann constant in eV/K (CODATA 2022)
constexpr double a0_m = 0.529177210544E-10; ///< Bohr radius in m (CODATA 2022)
constexpr double a0_cm = 0.529177210544E-8; ///< Bohr radius in cm (CODATA 2022)
constexpr double a0_nm = 0.0529177210544;   ///< Bohr radius in nm (CODATA 2022)
constexpr double eps0_As_Vm = 8.8541878188E-12; ///< vacuum electric permittivity in A s / V / m (CODATA 2022)
constexpr double au_F_V_cm = 5.1422082e9; ///< atomic unit of field strength in V/cm
constexpr double au_t_s = 2.4188843e-17; ///< atomic unit of time in s
} // namespace hydroconst
