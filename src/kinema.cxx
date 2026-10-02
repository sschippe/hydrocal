/**
 * @file kinema.cxx
 *
 * @brief Coding of relativistic kinematics of merged electron-ion beams
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: kinema.cxx 2037 2026-07-17 15:21:56Z iamp $
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */

#include "Cooler.h"
#include "hydroconst.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Relativistic formula for calculating the center-of-mass electron-ion
 * collision energy from the space-charge corrected lab electron energy
 *
 * @param Esc space-charge-corrected lab electron energy in eV
 * @param CoolingEnergy cooling energy in eV
 * @param IonMass ion mass in u
 * @param cosphi cosine of angle between electron beam and ion beam (default
 * value 1.0)
 *
 * @return center-of-mass electron-ion collision energy in eV, positive
 * (negative) if electrons are faster (slower) than ions
 */
double EcmFromEsc(double Esc, double CoolingEnergy, double IonMass,
                  double cosphi) {
  const double mec2 =
      hydroconst::mec2_eV; ///< electron rest mass times c2 in eV
  const double muc2 = hydroconst::muc2_eV; /// <atomic mass unit times c2 in eV
  double mic2 = IonMass * muc2;
  double mm = mec2 / mic2;
  double mm1 = 1.0 + mm;

  double gamma_e = 1.0 + Esc / mec2;
  double gamma_i = 1.0 + CoolingEnergy / mec2;
  double beta_e = sqrt(1.0 - 1.0 / (gamma_e * gamma_e));
  double beta_i = sqrt(1.0 - 1.0 / (gamma_i * gamma_i));
  double Gamma = gamma_i * gamma_e * (1.0 - beta_i * beta_e * cosphi);
  double sGamma = sqrt(1.0 + mm * (2.0 * Gamma + mm)) / mm1;
  double sign = Esc < CoolingEnergy ? -1.0 : 1.0;
  // return sign*mrc2*(Gamma-1.0);  // approximat Erel
  return sign * mic2 * mm1 * (sGamma - 1.0); // exact Ecm
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Relativistic formula for calculating the space-charge corrected lab
 * electron energy from the center-of-mass electron-ion collision energy
 *
 * @param Ecm center-of-mass electron-ion collision energy in eV
 * @param CoolingEnergy cooling energy in eV
 * @param IonMass ion mass in u
 * @param cosphi cosine of angle between electron beam and ion beam (default
 * value 1.0)
 *
 * @return space-charge-corrected lab electron energy in eV
 */
double EscFromEcm(double Ecm, double CoolingEnergy, double IonMass,
                  double cosphi) {
  const double mec2 =
      hydroconst::mec2_eV; ///< electron rest mass times c2 in eV
  const double muc2 = hydroconst::muc2_eV; /// <atomic mass unit times c2 in eV
  double mic2 = IonMass * muc2;
  // double mrc2 = mic2*mec2/(mic2+mec2);
  double mm = mec2 / mic2;
  double mm1 = 1.0 + mm;

  double sign = Ecm < 0.0 ? -1.0 : 1.0;
  double gamma_i = 1.0 + CoolingEnergy / mec2;
  double mmEcm = mm1 + fabs(Ecm) / mic2;
  double Gamma = 1.0 + 0.5 * (mmEcm * mmEcm - mm1 * mm1) / mm; // from exact Ecm
  double beta_i = sqrt(1.0 - 1.0 / (gamma_i * gamma_i));
  double one = gamma_i * gamma_i *
               (1.0 - beta_i * beta_i * cosphi * cosphi); // =1 for cosphi = 1
  double gamma_e =
      gamma_i * (Gamma + sign * beta_i * cosphi * sqrt(Gamma * Gamma - one)) /
      one;
  return mec2 * (double(gamma_e) - 1.0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Calculates space-charge corrected lab electron energy from the lab
 * electron energy without space-charge correction
 *
 * @param Elab lab electron energy in eV
 * @param Ie electron current in A
 * @param TubeBeamRatio vacuum-tube-radius divided by electron-beam-radius
 *
 * @return space-charge-corrected lab electron energy in eV
 */
double EscFromElab(double Elab, double Ie, double TubeBeamRatio) {
  const double Pi = hydroconst::pi; ///< Pi
  const double mec2 =
      hydroconst::mec2_eV; ///< electron rest mass times c2 in eV
  const double Clight = hydroconst::clight_m_s; ///< velocity of light in m/s
  const double eps0 =
      hydroconst::eps0_As_Vm; ///< permittivity of the vacuum in mA s/V m

  const unsigned int maxiter = 100;
  const double eps = 1.E-6;

  if ((Ie < 1E-6) || (TubeBeamRatio < 1.0))
    return Elab; // no space charge from zero electron current
  double factor = 0.25 * (1.0 + 2.0 * log(TubeBeamRatio)) /
                  (Pi * eps0 * Clight) * Ie / mec2;
  double Uk = Elab / mec2;
  double enew = 0.92 * Uk;
  double acc = 1.0;
  unsigned int i = 0;
  do {
    i++;
    double eold = enew;
    enew = Uk - factor * (eold + 1) / sqrt(eold * (eold + 2));
    acc = fabs(eold / enew - 1.0);
    if (i > maxiter) {
      cout << "Space-charge calculation 1:: Max number of iterations exceeded! "
              "Current accuracy: "
           << acc << endl;
      break;
    }
  } while (acc > eps);

  return enew * mec2;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Calculates the lab electron energy without space-charge correction
 * from the space-charge-corrected lab electron energy
 *
 * @param Esc space-charge-corrected lab electron energy in eV
 * @param Ie electron current in A
 * @param TubeBeamRatio vacuum-tube-radius divided by electron-beam-radius
 *
 * @return lab electron energy (without space-chare correction) in eV
 */
double ElabFromEsc(double Esc, double Ie, double TubeBeamRatio) {
  const unsigned int maxiter = 1000;
  const double eps = 1.E-6;

  if ((Ie < 1E-6) || (TubeBeamRatio < 1.0))
    return Esc;
  // the energy without space charge is higher than the energy with space charge
  double Edelta = EscFromElab(Esc, Ie, TubeBeamRatio) - Esc;
  double Elo = Esc - 2.0 * Edelta;
  double Ehi = Esc;
  double E_new = 0.5 * (Ehi + Elo);
  unsigned int i = 0;
  while (fabs(Ehi - Elo) > eps) {
    i++;
    E_new = 0.5 * (Ehi + Elo);
    double Esc_new = EscFromElab(E_new, Ie, TubeBeamRatio);
    // printf("E_new=%10.3f Esc_new=%10.3f Esc=%10.3f\n",E_new, Esc_new, Esc);
    if (Esc_new < Esc) {
      Ehi = E_new;
    } else {
      Elo = E_new;
    }
    if (i > maxiter) {
      cout << "Space-charge calculation 2:: Max number of iterations exceeded! "
              "Current accuracy: "
           << fabs(Ehi - Elo) << endl;
      break;
    }
  }
  return E_new;
}

///////////////////////////////////////////////////////////////////////////

void kinema(void) {
  const double mec2 = hydroconst::mec2_eV;
  const double muc2 = hydroconst::muc2_eV;
  ;

  cout << endl
       << endl
       << "**** Kinematics in electron-ion merged-beams and crossed-beams "
          "experiments ****"
       << endl;
  int storage_ring_id =
      COOLER::SelectStorageRing(true);  // see enum StorageRingIDs in Cooler.h
  COOLER cooler(storage_ring_id, true); // initializes cooler dimensions
  double Ucool = cooler.CoolingVoltage();
  double Ie = cooler.ElectronCurrent();
  double tube_beam_ratio = cooler.TubeBeamRatio();

  double Ecool = EscFromElab(Ucool, Ie, tube_beam_ratio);
  double gamma_i = 1.0 + Ecool / mec2;
  double beta_i = sqrt(1.0 - 1.0 / (gamma_i * gamma_i));
  double Ei = Ecool / mec2 * muc2 * 1E-6; // in MeV/u
  double Doppler_factor = sqrt((1.0 + beta_i) / (1.0 - beta_i));

  cout << endl;
  cout << " cooling voltage (V) : " << Ucool << endl;
  cout << " cooling energy (eV) : " << Ecool << endl;
  cout << " space charge (eV) ..: " << Ucool - Ecool << endl;
  cout << " ion energy (MeV/u) .: " << Ei << endl;
  cout << " ion gamma ..........: " << gamma_i << endl;
  cout << " ion beta ...........: " << beta_i << endl;
  cout << " Doppler factor .....: " << Doppler_factor << endl;

  double A, Udmin, Udmax, Udelta, costheta;
  cout << endl << " Give ion mass in u .....................................: ";
  cin >> A;
  if (storage_ring_id == COOLER::CROSSED_BEAMS) {
    cout << endl
         << " Give min and max e-gun cathode voltage in V ............: ";
    costheta = 0.0;
  } else {
    cout << endl
         << " Give min and max detuning voltage in V .................: ";
    costheta = 1.0;
  }
  cin >> Udmin >> Udmax;
  cout << endl << " Give step with in eV ...................................: ";
  cin >> Udelta;

  printf("\n    Ud (V)     Ue (V)   Esc (eV)    Ecm (eV)\n");
  for (double Ud = Udmin; Ud <= Udmax; Ud += Udelta) {
    double Ue = Ud;
    if (storage_ring_id != COOLER::CROSSED_BEAMS)
      Ue += Ucool;
    double Esc = EscFromElab(Ue, Ie, tube_beam_ratio);
    double Ecm = EcmFromEsc(Esc, Ecool, A, costheta);
    printf("%10.3f %10.3f %10.3f  %10.5f\n", Ud, Ue, Esc, Ecm);
  }
}

///////////////////////////////////////////////////////////////////////////

double weizmass(int A, int Z, int Q) {
  // returns nuclear mass for Q>=0 or
  // binding energy per nucleon for Q<0

  if ((Z < 1) || (A < 1) || (Z > A))
    return 0.0;

  const double aV = 15.85e6;
  const double aS = 18.34e6;
  const double aC = 0.71e6;
  const double aA = 92.86e6;
  const double aP = 11.46e6;
  const double mp = hydroconst::mpc2_eV;
  const double mn = hydroconst::mnc2_eV;
  const double me = hydroconst::mec2_eV;

  double mH = mp + me - hydroconst::Ryd_eV;
  double A3, B, Tz;
  int N, sign, Zeven, Neven;

  N = A - Z;
  Tz = 0.5 * (Z - N);

  Zeven = Z % 2;
  Neven = N % 2;
  sign = (Zeven == Neven) ? (Zeven ? 1 : -1) : 0;

  A3 = exp(log(double(A)) / 3.0);
  B = aV * A - aS * A3 * A3 - aC * Z * Z / A3 - aA * Tz * Tz / A +
      sign * aP / sqrt(double(A));

  double wm = (Q >= 0) ? Z * mH + N * mn - Q * me - B : B / double(A);

  return wm;
}

///////////////////////////////////////////////////////////////////////////

void TestWeizmass(void) {
  int a, z, q, mode, amax, amin;
  double wm;
  ofstream fout;
  string fn;

  const double malpha = hydroconst::malphac2_eV;
  printf("\n Give maximum mass number : ");
  scanf("%d", &amax);
  printf("\n Mode 0 : max(B(Z,A)/A) as function of A");
  printf("\n Mode 1 : B(Z,A)/A");
  printf("\n Mode 2 : m(Z,A)");
  printf("\n Mode 3 : m(Z,A)-m(Z-2,A-4)-m(alpha)");
  printf("\n Give mode ...............: ");
  scanf("%d", &mode);
  printf("\n Give filename for output : ");
  cin >> fn;

  fout.open(fn);
  amin = (mode < 3) ? 1 : 5;
  q = (mode > 1) ? 0 : -1;
  for (a = amin; a <= amax; a++) {
    int zmax = 1;
    double wmax = 0.0;
    for (z = amin; z <= amax; z++) {
      wm = weizmass(a, z, q);
      if (wm > wmax) {
        wmax = wm, zmax = z;
      }
      if (mode > 2) {
        double wm2 = weizmass(a - 4, z - 2, q);
        if (wm2 > 0) {
          wm -= wm2 + malpha;
        } else {
          wm = -malpha;
        }
      }
      if (mode > 0)
        fout << setw(12) << setprecision(6) << wm / 1e6 << " ";
    }
    if (mode > 0) {
      fout << "\n";
    } else {
      double wzmax = weizmass(a, zmax, 0);
      double ealpha = wzmax - weizmass(a - 4, zmax - 2, 0) - malpha;
      double ecoul = -1.11e6 * z * z * exp(-log(0.5 * a) / 3) / 8.0;
      double efiss = weizmass(a / 2, zmax / 2, 0);
      efiss = (efiss > 0) ? wzmax - 2 * efiss : 0.0;
      fout << setw(5) << a << " " << setw(5) << zmax << " " << setw(12)
           << setprecision(6) << wmax / 1e6 << " " << setw(12)
           << setprecision(6) << wzmax / 1e6 << " " << setw(12)
           << setprecision(6) << ealpha / 1e6 << " " << setw(12)
           << setprecision(6) << efiss / 1e6 << " " << setw(12)
           << setprecision(6) << ecoul / 1e6 << "\n";
    }
  }
  fout.close();
  cout << "\n Mass output  written to file " << fn << ".\n\n";
}

///////////////////////////////////////////////////////////////////////////
