/**
 * @file RRDRconvolve.cxx
 *
 * @brief Recombination rate coefficients from tabulated RR or DR cross sections
 *
// SPDX-License-Identifier: MIT
 *
 * The convolution is with
 * - either an isotropic Maxwellian yield a plasma rate coefficient as function
of electron temperature [convolute_alpha()]
 * - or a flattened or gaussian energy distribution yielding a convoluted cross
section and a rate  coefficient [convflat()].
 */

#include "RRDRxsec.h"
#include "fele.h"
#include "hydroconst.h"
#include "hydromath.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "stdin_guard.h"

using namespace std;

/**
 * @brief Convolution yielding a plasma rate coefficient
 *
 * The convolution is carried out numerically. Input data do not have to be
 * spaced equidistantly.
 */
void RRDRconvolve_PlasmaRateCoeff(void) {
  string line, filename;
  ifstream fin, ftemp;
  ofstream fout;
  int n = 0, choice = 0;
  double elo = 0.0, ehi = 0.0;
  cout << "\n Convolution of cross sections with a Maxwellian electron energy";
  cout << "\n distribution yielding the rate coefficient in a plasma as a";
  cout << "\n function of electron temperature kT (in eV).\n";

  cout << "\n Cross section data are read in from files containing at least";
  cout << "\n two columns. The first one is assumed to contain energies";
  cout << "\n in eV and the second either cross sections in cm^2 or rate";
  cout << "\n coefficients in cm^3/s. Before being convoluted the latter";
  cout << "\n are converted in to cross sections by division by the electron";
  cout << "\n velocity. The convolution is carried out numerically. Energies";
  cout << "\n do not have to be spaced equidistantly.\n";

  cout << "\n The calculation of semiclassical RR rate coefficients with";
  cout << "\n Stobbe corrections from the corresponding cross section";
  cout << "\n is carried out analytically without numerical convolution.\n";

  cout << "\n Convolute what?";
  cout << "\n   1: rate coefficients from file";
  cout << "\n   2: cross sections from file (electron ion collisions)";
  cout << "\n   3: semiclassical RR cross section with Stobbe corrections";
  cout << "\n   4: cross sections from file (nulear reactions) \n";
  while ((choice < 1) || (choice > 4)) {
    cout << "\n Make a choice .............................................: ";
    cin >> choice;
  }

  double m1 = hydroconst::mec2_eV; // projectile mass
  double mm = 1.0;                 // mass ratio (reduced mass)/m1;
  int xsec_mode = 0;
  switch (choice) {
  case 4:
    xsec_mode = 1;
    double m2;
    cout << "\n Give target mass in atomic mass units ....: ";
    cin >> m2;
    cout << "\n Give projectile mass in atomic mass units : ";
    cin >> m1;
    mm = m2 / (m1 + m2);       // mass ratio (reduced mass)/m1;
    m1 *= hydroconst::muc2_eV; // conversion to eV
    cout << "\n Give name of file containing";
    cout << " energies (in eV) and cross sections ..: ";
    cin >> filename;
    break;
  case 3:
    calcAlphaRRplasma();
    return;
    break;
  case 2:
    xsec_mode = 1;
    cout << "\n Give name of file containing";
    cout << " energies (in eV) and cross sections ..: ";
    cin >> filename;
    break;
  case 1:
    cout << "\n Give name of file containing";
    cout << " energies (in eV) and rate coefficients: ";
    cin >> filename;
    break;
  }

  fin.open(filename);
  while (getline(fin, line))
    n++;
  fin.close();
  int nlines = n;
  std::vector<double> ecm(nlines);
  std::vector<double> alpha(nlines);
  fin.open(filename);
  for (n = 0; n < nlines; n++) {
    getline(fin, line);
    istringstream ss(line);
    ss >> ecm[n] >> alpha[n];
    cout << setw(10) << uppercase << defaultfloat << setprecision(4) << ecm[n]
         << " " << setw(10) << uppercase << defaultfloat << setprecision(4)
         << alpha[n] << "\n";
  }
  fin.close();
  cout << "\n " << setw(5) << nlines << " lines read !\n";

  if (ehi <= elo) {
    elo = ecm[0];
    ehi = ecm[nlines - 1];
  }
  double tmin, tmax, dt, t, kt;
  cout << "\n Give name of output file ..................................: ";
  cin >> filename;
  fout.open(filename);
  fout << "         T           kT         alpha alpha*kt^3/2\n\n";

  int steps;
  cout << "\n Read temperatures from file? Give filename (0=no file) ....: ";
  cin >> filename;
  int read_temp = (filename != "0");
  if (read_temp) {
    ftemp.open(filename);
    if (!ftemp.is_open()) {
      read_temp = 0;
      cout << "\n file " << filename << " not found!";
    } else { // count the number of energies given in the data file
      steps = 0;
      while (getline(ftemp, line))
        steps++;
      steps--;
      ftemp.close();
    }
  }

  if (read_temp) {
    ftemp.open(filename);
  } else {
    cout << "\n Give log min, max, delta kT in eV (separated by spaces) ...: ";
    cin >> tmin >> tmax >> dt;
    steps = int((tmax - tmin) / dt);
  }

  for (int i = 0; i <= steps; i++) {
    if (read_temp) {
      getline(ftemp, line);
      istringstream ss(line);
      ss >> kt;
    } else {
      t = tmin + i * dt;
      kt = exp(log(10.0) * t);
    }
    double delta = 0.5 * (ecm[1] - ecm[0]);
    double f = sqrt(ecm[0]) * exp(-ecm[0] * mm / kt);
    if (xsec_mode)
      f *= sqrt(ecm[0]);
    double sum = ecm[0] < elo ? 0.0 : delta * f * alpha[0];
    for (n = 1; n < nlines - 1; n++) {
      if ((elo < ecm[n]) && (ecm[n] < ehi)) {
        delta = 0.5 * (ecm[n + 1] - ecm[n - 1]);
        f = sqrt(ecm[n]) * exp(-ecm[n] * mm / kt);
        if (xsec_mode)
          f *= sqrt(ecm[n]);
        sum += delta * f * alpha[n];
      }
    }
    if (ecm[nlines - 1] <= ehi) {
      delta = 0.5 * (ecm[nlines - 1] - ecm[nlines - 2]);
      f = sqrt(ecm[nlines - 1]) * exp(-ecm[nlines - 1] * mm / kt);
      if (xsec_mode)
        f *= sqrt(ecm[nlines - 1]);
      sum += delta * f * alpha[nlines - 1];
    }
    sum *= 2.0 * mm * sqrt(mm / hydroconst::pi);
    if (xsec_mode)
      sum *= sqrt(2.0 / m1) * hydroconst::clight_cm_s * 100.0;
    fout << setw(10) << uppercase << defaultfloat << setprecision(4)
         << kt / hydroconst::kB_eV_K << "   " << setw(10) << uppercase
         << defaultfloat << setprecision(4) << kt << "    " << setw(10)
         << uppercase << defaultfloat << setprecision(4) << sum / kt / sqrt(kt)
         << "   " << setw(10) << uppercase << defaultfloat << setprecision(4)
         << sum << "\n";
  }
  fout.close();
  if (read_temp)
    ftemp.close();
}

/***************************************************************************/
/**
 * @brief Convolution yielding a convolved cross section and a merged-beams rate
 * coefficient
 *
 * @param data_type = 1: convolve cross section data
 * @param data_type = 2: convolve rate coefficient data
 *
 * The convolution is carried out numerically. Input data do not have to be
 * spaced equidistantly.
 */
void RRDRconvolve_MBrateCoeff(int datatype) {
  if (datatype == 1) {
    cout << "\n Convolution of cross sections contained as x-, y- columns";
    cout << "\n separated by a <space> in a file ";
  }
  if (datatype == 2) {
    cout << "\n Derivation of an experimental cross section from a rate "
            "coefficient";
    cout << "\n contained as <space> separated  x-, y- columns in a file and "
            "convolution";
    cout << "\n of the cross sections ";
  }
  cout << "with a flattened Maxwellian, with a Gaussian, or a";
  cout << "\n normalized trapezoid electron energy distribution yielding";
  cout << "\n a convoluted cross section and a merged beams rate coefficient.";
  cout << "\n \n";

  if (datatype == 1) {
    cout << "\n Cross section data are read in from files containing at least";
  }
  if (datatype == 2) {
    cout << "\n Rate coefficient data are read in from files containing at "
            "least";
  }
  cout << "\n two columns. The first one is assumed to contain energies";
  if (datatype == 1) {
    cout << "\n and the second cross sections.\n";
  }
  if (datatype == 2) {
    cout << "\n and the second rate coefficients.\n";
  }

  cout << "\n The convolution is carried out numerically. Energies";
  cout << "\n do not have to be spaced equidistantly.\n";

  cout << "\n The output file contains four columns: ";
  cout << "\n Energy (eV), convoluted cross section (cm^2),";
  cout << "\n rate coeffient (cm^3/s), and convolution width (eV).\n";

  int choice;
  cout << "\n Convolute with?";
  cout << "\n   1: flattened Maxwellian (or isotropic with FWHM ~ E^0.5)";
  cout << "\n   2: normalized Gaussian (FWHM energy independent)";
  cout << "\n   3: normalized trapezoid";
  cout << "\n      Make a choice ......................................: ";
  cin >> choice;

  string filename, filenameout, line;
  ifstream fin, fenergy;
  ofstream fout;
  double e_unit, cs_unit, tolerance;
  int fin_exists = 0;

  while (!fin_exists) {
    cout << "\n Unconvoluted cross sections are now read from a file.";
    cout << "\n ATTENTION: The file must not contain empty lines!";
    cout << "\n Give name of input file..................................: ";
    cin >> filename;
    fin.open(filename);
    if (fin.is_open())
      fin_exists = 1;
  }

  int nlines = 0, n;
  while (getline(fin, line)) {
    if (line[0] == '#') {
      cout << "# " << line << "\n";
    } else {
      nlines++;
    }
  }
  fin.close();

  cout << "\n Give energy unit in eV ..................................: ";
  cin >> e_unit;
  cout << "\n Give cross section unit in cm^2 .........................: ";
  cin >> cs_unit;

  cout << "\n\n If the x-axis spacing is too coarse the convolution is bound "
          "to fail.";
  cout << "\n In this case the convoluted cross section is substitued by the "
          "original.";
  cout << "\n cross section. The associated quality of the convolution is "
          "monitored by";
  cout << "\n checking whether the normalization integral of the distribution";
  cout << "\n function evaluates to 1. The user can allow for tiny deviations "
          "from 1 by";
  cout << "\n making an appropriate choice for the 'tolerance'. A reasonable "
          "value would";
  cout << "\n be 1e-4. The check can be switched off by setting 'tolerance <= "
          "0.0'.\n";

  cout << "\n Give tolerance (<=0 means no check) .....................: ";
  cin >> tolerance;

  double elo = 1E99, ehi = -1E99;
  vector<double> ecm(nlines, 0.0);
  vector<double> xsec(nlines, 0.0);

  fin.open(filename);
  for (n = 0; n < nlines; n++) {
    getline(fin, line);
    if (line[0] == '#') {
      n--;
      continue;
    } // skip comment lines
    istringstream ss(line);
    ss >> ecm[n] >> xsec[n];
    if (ecm[n] < 0) {
      nlines--;
      n--;
      continue;
    }
    ecm[n] *= e_unit; // energies are now in eV
    if (datatype == 2) {
      xsec[n] /=
          (hydroconst::clight_cm_s * sqrt(2.0 * ecm[n] / hydroconst::mec2_eV));
    }
    xsec[n] *= cs_unit; // cross sections are now in cm^2
    if (ecm[n] < elo)
      elo = ecm[n];
    if (ecm[n] > ehi)
      ehi = ecm[n];
  }
  fin.close();
  cout << "\n Unconvoluted cross sections read at " << nlines + 1
       << " energies";
  cout << "\n                             ranging from " << uppercase
       << defaultfloat << setprecision(6) << elo << " to " << uppercase
       << defaultfloat << setprecision(6) << ehi << " eV.\n";

  // calculate second derivative for spline interpolation below
  vector<double> dderiv(nlines, 0.0);
  spline(nlines, ecm, xsec, dderiv);

  FELE fele = NULL;
  double fwhm = 0.0, ktpar, ktperp;
  switch (choice) {
  case 1:
    cout << "\n Give parallel electron beam temperature in meV ..........: ";
    cin >> ktpar;
    ktpar *= 0.001;
    cout << "\n Give perpendicular electron beam temperature in meV .....: ";
    cin >> ktperp;
    ktperp *= 0.001;
    fele = &fecool;
    break;
  case 2:
    cout << "\n Give fwhm in eV .........................................: ";
    cin >> ktpar;
    fele = &felegauss;
    break;
  case 3:
    cout << "\n Give base width in eV ...................................: ";
    cin >> ktpar;
    cout << "\n  Give top width in eV ...................................: ";
    cin >> ktperp;
    fele = &trapezoid;
    break;
  }

  double eV, emin, emax, edelta;
  int steps;
  // in order to generate theoretical rate coefficients at exactly the
  // same energies as given in an experimental data file, the energies
  // from this file can be read in. It is assumed that energies are
  // listed in the first column of the experimental data file.
  cout << "\n Read energies from file? Give filename (0=no file) ......: ";
  line = "0";
  cin >> filename;
  int read_energy = (filename != "0");
  if (read_energy) {
    fenergy.open(filename);
    if (!fenergy.is_open()) {
      read_energy = 0;
      cout << "\n file " << filename << " not found!";
    } else { // count the number of energies given in the data file
      steps = 0;
      while (getline(fenergy, line))
        steps++;
      steps--;
      fenergy.close();
    }
  }
  if (!read_energy) {
    cout << "\n Give energy range (emin,emax,delta) .....................: ";
    cin >> emin >> emax >> edelta;
    steps = int((emax - emin) / edelta);
  }
  cout << "\n Convolved cross sections will be calculated at " << steps
       << " energies.\n";
  cout << "\n Give name of output file ................................: ";
  cin >> filenameout;
  cout << "\n Progress: ";
  fout.open(filenameout);
  fout << "      Erel           sigma           alpha                  FWHM\n";
  fout << "      (eV)           (cm^2)          (cm^3/s)               (eV)\n";

  if (read_energy) {
    fenergy.open(filename);
    if (!fenergy.is_open())
      cerr << "Error opening file";
  }

  for (int i = 0; i <= steps; i++) { // loop over energies
    if ((i % 100) == 0)
      cout << "." << flush;
    if (read_energy) {
      getline(fenergy, line);
      istringstream ss(line);
      ss >> eV;
    } else {
      eV = emin + i * edelta;
      if (emin < 0)
        eV = exp(log(10.0) * eV);
    }

    double fpar, fperp;
    switch (choice) {
    case 1:
      fpar = 4 * sqrt(ktpar * eV * log(2.0));
      fperp = log(2.0) * ktperp;
      fwhm = sqrt(fpar * fpar + fperp * fperp);
      break;
    case 2:
      fwhm = ktpar;
      break;
    case 3:
      fwhm = 0.5 * (ktpar + ktperp);
      break;
    }

    double delta = ecm[1] - ecm[0];
    double df = delta * fele(eV, ecm[0], ktpar, ktperp);
    double sumdf = ecm[0] < elo ? 0.0 : df;
    double sumdfx = ecm[0] < elo ? 0.0 : df * xsec[0];

    for (n = 1; n < nlines - 1; n++) {
      if ((elo < ecm[n]) && (ecm[n] < ehi)) {
        delta = 0.5 * (ecm[n + 1] - ecm[n - 1]);
        df = delta * fele(eV, ecm[n], ktpar, ktperp);
        sumdf += df;
        sumdfx += df * xsec[n];
      }
    }

    if (ecm[nlines - 1] <= ehi) {
      delta = ecm[nlines - 1] - ecm[nlines - 2];
      df = delta * fele(eV, ecm[n], ktpar, ktperp);
      sumdf += df;
      sumdfx += df * xsec[nlines - 1];
    }

    // Here we check wether the integration yields a normalized distribution
    // function. This is not the case if the cross-section grid is too coarse,
    // i.e. when the width of the distribution function is smaller than the
    // difference between adjacent energy points of the cross-section grid. In
    // this case we substitute the convoluted cross section by the cross section
    // itself.
    if ((tolerance > 0.0) && (fabs(sumdf - 1.0) > tolerance)) {
      // cubic spline interpolation

      sumdfx = splint(eV, nlines, ecm, xsec, dderiv);
    }

    if (fabs(sumdfx) < 1E-99)
      sumdfx = 0.0;
    double rate_coeff =
        sumdfx * sqrt(2 * eV / hydroconst::mec2_eV) * hydroconst::clight_cm_s;
    fout << setw(15) << uppercase << defaultfloat << setprecision(8) << eV
         << "   " << setw(15) << uppercase << defaultfloat << setprecision(8)
         << sumdfx << "   " << setw(15) << uppercase << defaultfloat
         << setprecision(8) << rate_coeff << "   " << setw(15) << uppercase
         << defaultfloat << setprecision(8) << fwhm << "\n";
  } // end for(i...)
  fout.close();
  if (read_energy)
    fenergy.close();
  cout << "\n\n";
}
