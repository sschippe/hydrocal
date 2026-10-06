// SPDX-License-Identifier: MIT
/**
 * @file peakfunctions.h
 * @brief Definitions of peak functions
 *
 * @author Stefan Schippers
 * @date 2013-04-23
 *
 */

#pragma once

enum PEAKSHAPES {
  PEAK_undefined,
  PEAK_delta,
  PEAK_Lorentz,
  PEAK_Fano,
  PEAK_Gauss,
  PEAK_Voigt,
  PEAK_FanoVoigt,
  PEAK_Lorentz_Steih,
  PEAK_Lorentz_Amaro,
  PEAK_P04,
  PEAK_P04_Lorentz
};

void sum_of_peaks(void);

double voigt(double E, double E0, double a, double wl, double wg);

double gauss(double E, double E0, double a, double wg);

double lorentz(double E, double E0, double a, double wl);

double fano(double E, double E0, double a, double q, double wl);

double fanovoigt(double E, double E0, double a, double q, double wl, double wg);

double fanovoigtderiv(double E, double E0, double a, double q, double wl,
                      double wg, double &df_dE0, double &df_da, double &df_dq,
                      double &df_dwl, double &df_dwg);

double trapezlorentzian(double E, double E0, double a, double wl, double wb,
                        double wt);

double RydbergSeries(double x, double dx, double wG, double zeff, double slim,
                     double m0, double m1, double q0, double q1, double w0,
                     double w1);

double P04ped(double E, double Enom, double exitslit);

double P04lorentz(double E, double Eres, double S, double wL, double exitslit,
                  double &error, double &Eintmin, double &Eintmax);
