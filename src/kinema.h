/**
 * @file kinema.h
 *
 * @brief Exports functions related to relativistic kinematics of merged
electron-ion beams
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */

#pragma once

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Relativistic formula for calculating the space-charge corrected lab
 * electron energy from the center-of-mass electron-ion collision  energy
 */
double EscFromEcm(double Ecm, double CoolingEnergy, double IonMass,
                  double cosphi = 1.0);

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Relativistic formula for calculating the center-of-mass electron-ion
 * collision energy from the space-charge correted lab electron energy
 */
double EcmFromEsc(double Esc, double CoolingEnergy, double IonMass,
                  double cosphi = 1.0);

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Calculates space-charge corrected lab electron energy from the lab
 * electron energy without space-charge correction
 */
double EscFromElab(double Elab, double Ie, double TubeBeamRatio);

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Calculates the lab electron energy without space-charge correction
 * from the space-charge-corrected lab electron energy
 */
double ElabFromEsc(double Esc, double Ie, double TubeBeamRatio);

void kinema(void);

void TestWeizmass(void);
