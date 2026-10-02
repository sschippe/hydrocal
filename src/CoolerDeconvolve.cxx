/**
 * @file CoolerDeconvolve.cxx
 *
 * @brief Toroid deconvolution of experimental merged-beams rate coeffcients
 *
 * Implements algorithm of Lampert et al., Phys. Rev. A 53, 1413–1423 (1996);
 * https://doi.org/10.1103/PhysRevA.53.1413
 *
 * @author Stefan Schippers
 * @verbatim
 $Id: CoolerDeconvolve.cxx 2037 2026-07-17 15:21:56Z iamp $
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
#include <chrono> // clocks and time
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief center-of-mass in the cooler from nominal center-of-mass energy
 *
 * @param Ecm on entry: nominal center-of-mass energy, on exit center-of-mass
 * energy accounting for deflection angle in cooler
 * @param z position in the cooler in cm (z=0 in the center of the cooler)
 * @param cooler the electron cooler in use
 * @param cooling_energy cooling energy in eV
 * @param ion_mass ion mass in u
 *
 * @return true if there is overlap with the ion beam
 * @return false if there is no overlap with the ion beam
 */
bool EcmCooler(double &Ecm, double z, COOLER cooler, double cooling_energy,
               double ion_mass) {
  double theta, y;
  cooler.BeamAngle(z, theta, y);
  if (y > cooler.BeamRadius()) {
    return false; // idicates that there is no overlap
  } else {
    double Esc = EscFromEcm(Ecm, cooling_energy, ion_mass, 1.0);
    Ecm = EcmFromEsc(Esc, cooling_energy, ion_mass, cos(theta));
    return true;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Toroid deconvolution of experimental merged-beams rate coefficients
 *
 * see A. Lampert et al., Phys. Rev. A 53, 1413 (1996)
 * (doi: 10.1103/PhysRevA.53.1413)
 */
void CoolerDeconvolve(void) {
  // *****************************  time stamp
  // **************************************
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  // ****************************** read data
  // *******************************************
  string filename;
  cout << endl
       << " **** Cooler deconvolution of experimental merged-beams rate "
          "coefficients ****"
       << endl;
  int storage_ring_id = COOLER::SelectStorageRing(); // see Cooler.h
  COOLER cooler(storage_ring_id, false);             // see Cooler.h
  double max_overlap_length = cooler.MaxOverlapLength();
  double nominal_overlap_length = cooler.NominalOverlapLength();
  double cooling_energy = cooler.CoolingEnergy();
  cout << endl;
  cout << endl
       << " The deconvolution requires an input file with experimental data.";
  cout << endl
       << " Each line must contain (separated by spaces) <Ecm (eV)> <alpha "
          "(cm3/s)> <error (cm3/s)>.";
  cout << endl
       << " The data are assumed to have been sorted in ascending order of the "
          "energies.";
  cout << endl << " Give the filename of input file ........................: ";
  cin >> filename;
  fstream fin;
  vector<double> ecm, alpha_0, error_0;
  double dEcm_min = 1E9;
  int npts = 0;
  fin.open(filename, fstream::in);
  while (!fin.fail() && !fin.eof()) {
    double val1, val2, val3;
    if (!(fin >> val1 >> val2 >> val3))
      break; // stop on failed extraction (e.g. trailing newline / EOF)
    ecm.push_back(val1);
    alpha_0.push_back(val2);
    error_0.push_back(val3);
    if (npts > 0) {
      if (ecm[npts] <= ecm[npts - 1]) {
        cout << endl
             << "ERROR: Energies are not in ascending order: " << ecm[npts - 1]
             << "  " << ecm[npts] << endl
             << endl;
        exit(0);
      }
      double dEcm = ecm[npts] - ecm[npts - 1];
      if (dEcm < dEcm_min)
        dEcm_min = dEcm;
    }
    npts++;
  }
  fin.close();
  cout << endl
       << " " << npts << " data points read from file " << filename << "."
       << endl;
  cout << endl
       << " The nominal overlap length should be the same as the one used for "
          "generating these data.";
  cout << endl
       << " Currently it is set to a value of " << nominal_overlap_length
       << " cm.";
  cout << endl << " Do you need to change this value? (y/n) ................: ";
  char answer;
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    cout << endl
         << " Give new value for nominal overlap length (cm) .........: ";
    cin >> nominal_overlap_length;
    cooler.SetNominalOverlapLength(nominal_overlap_length);
  }
  double ion_mass;
  cout << endl << " Give ion mass in u .....................................: ";
  cin >> ion_mass;
  cout << endl;
  double dz, eps;
  int maxiter;
  cout << endl
       << " The integration length along the cooler axis in maximally "
       << max_overlap_length << " cm.";
  cout << endl << " Give step size in electron cooler (cm) .................: ";
  cin >> dz;
  cout << endl << " Give relative accuracy goal ............................: ";
  cin >> eps;
  cout << endl << " Give maximum number of iterations.......................: ";
  cin >> maxiter;
  cout << endl << endl;

  // ****************************** setup iteration
  // **************************************
  vector<double> alpha_old, error_old, alpha_new, error_new;
  alpha_old.resize(npts);
  error_old.resize(npts);
  alpha_new.resize(npts);
  error_new.resize(npts);

  for (int n = 0; n < npts; n++) {
    error_0[n] *= error_0[n]; // errors are added in quadrature
    alpha_new[n] = alpha_0[n];
    error_new[n] = error_0[n];
  }
  int iteration_counter = 0;
  double max_rel_diff = 1;
  cout << " current iteration:";

  // ****************************** perform iteration
  // **************************************
  while (max_rel_diff > eps) {
    iteration_counter++;
    max_rel_diff = 0.0;
    alpha_old = alpha_new;
    error_old = error_new;
    for (int n = 0; n < npts; n++) {
      double alpha_cooler = 0.0, error_cooler = 0.0, overlap_length = 0.0;
      for (double zz = 0; zz <= 0.5 * max_overlap_length; zz += dz) {
        // integrate along cooler axis
        int nz = 0;
        double ecm_cooler = ecm[n];
        if (EcmCooler(ecm_cooler, zz, cooler, cooling_energy, ion_mass)) {
          overlap_length += dz;
          if (fabs(ecm_cooler - ecm[n]) < 0.1 * dEcm_min) {
            // energy did not change significantly
            nz = n;
            alpha_cooler += alpha_old[nz];
            error_cooler += error_old[nz];
          } else {
            // search for index of new energy
            for (nz = 0; nz < npts; nz++)
              if (ecm[nz] > ecm_cooler)
                break;
            if (nz == npts)
              nz = npts - 1; // if break has not occurred
            if ((nz == 0) || (nz == (npts - 1))) {
              alpha_cooler += alpha_old[nz];
              error_cooler += error_old[nz];
            } else {
              // interpolate alpha(ecm_z)
              double r = (ecm_cooler - ecm[nz - 1]) / (ecm[nz] - ecm[nz - 1]);
              alpha_cooler += alpha_old[nz - 1] * (1.0 - r) + r * alpha_old[nz];
              double error_interp =
                  sqrt(error_old[nz - 1]) * (1.0 - r) + r * sqrt(error_old[nz]);
              error_cooler +=
                  error_interp *
                  error_interp; // note that errors are added in quadrature
            }
            // cout << iteration_counter << "  "  << n << "   " << ecm[n] << "
            // " << ecm_cooler << "  " << nz << "  " << alpha_old[nz] << "  " <<
            // alpha_cooler << endl;
          } // end else if ( fabs(ecm_cooler ...)...)
        } // end if (ECMCooler(...))
      } // end for (zz...)
      double length_ratio1 = dz / overlap_length;
      alpha_cooler *= length_ratio1;
      error_cooler *= length_ratio1 * length_ratio1;
      double length_ratio2 = 0.5 * nominal_overlap_length / overlap_length;
      alpha_new[n] = alpha_old[n] - alpha_cooler + alpha_0[n] * length_ratio2;
      error_new[n] = error_old[n] + error_cooler +
                     error_0[n] * length_ratio2 * length_ratio2;
      double rel_diff = fabs(alpha_old[n]) > 0.0
                            ? fabs(1.0 - fabs(alpha_new[n] / alpha_old[n]))
                            : 0.0;
      if (rel_diff > max_rel_diff)
        max_rel_diff = rel_diff; // disregard end points
    } // end for (int n...)
    cout << "    " << iteration_counter << " (" << max_rel_diff << ")";
    cout.flush();
    if (iteration_counter >= maxiter) {
      cout << endl
           << endl
           << " ATTENTION: Max number of iterations (" << maxiter
           << ") exceeded. No convergence achieved." << endl;
      break;
    }
  } // end while (max_rel_dif>eps)

  // check that the convolution of the deconvolved rate coefficient yields again
  // the input rate coefficient (the errors will be wrong!)
  for (int n = 0; n < npts; n++) {
    double alpha_cooler = 0.0;
    for (double zz = 0; zz <= 0.5 * max_overlap_length; zz += dz) {
      // integrate along cooler axis
      int nz = 0;
      double ecm_cooler = ecm[n];
      if (EcmCooler(ecm_cooler, zz, cooler, cooling_energy, ion_mass)) {
        if (fabs(ecm_cooler - ecm[n]) < 0.1 * dEcm_min) {
          // energy did not change significantly
          nz = n;
          alpha_cooler += alpha_new[nz];
        } else {
          // search for index of new energy
          for (nz = 0; nz < npts; nz++)
            if (ecm[nz] > ecm_cooler)
              break;
          if (nz == npts)
            nz = npts - 1; // if break has not occurred
          if ((nz == 0) || (nz == (npts - 1))) {
            alpha_cooler += alpha_new[nz];
          } else {
            // interpolate alpha(ecm_z)
            double r = (ecm_cooler - ecm[nz - 1]) / (ecm[nz] - ecm[nz - 1]);
            alpha_cooler += alpha_new[nz - 1] * (1.0 - r) + r * alpha_new[nz];
          }
        } // end else if ( fabs(ecm_cooler ...)...)
      } // end if (ECMCooler(...))
    } // end for (zz...)
    double length_ratio =
        2.0 * dz / nominal_overlap_length; // factor 2 because we step only
                                           // through one half of the cooler
    alpha_old[n] = alpha_cooler * length_ratio;
    error_old[n] = (alpha_0[n] > 0.0)
                       ? fabs(1.0 - alpha_old[n] / alpha_0[n])
                       : -1.0; // relative deviation between test and original
  } // end for (int n...)

  // ****************************** output of results
  // **************************************

  fstream fout;
  string outfilename = filename;
  size_t dotpos = outfilename.rfind('.');
  if (dotpos != string::npos)
    outfilename.erase(dotpos + 1,
                      outfilename.length() - dotpos); // truncate extension
  else
    outfilename.append("."); // no extension present: append separator
  outfilename.append("dcv");
  fout.open(outfilename, fstream::out);
  fout << "####################################################################"
          "##########################"
       << endl;
  fout << "### Cooler deconvolution of rate-coefficient data from an "
          "electron-ion merged-beams experiment"
       << endl;
  fout << "###" << endl;
  fout << "###  hydrocal revision     : " << HYDROCAL_REVISION
       << endl; // defined in buildinfo.h
  fout << "###               filename : " << outfilename << endl;
  fout << "###         input filename : " << filename << endl;
  fout << "###    number of data sets : " << npts << endl;
  fout << "###      start date & time : " << asctime(timeinfo);
  cooler.Print(fout);
  fout << "###    cooling energy (eV) : " << cooling_energy << endl;
  fout << "###           ion mass (u) : " << ion_mass << endl;
  fout << "### min Ecm step size (cm) : " << dEcm_min << endl;
  fout << "### integr. step size (cm) : " << dz << endl;
  fout << "###          accuracy goal : " << eps << endl;
  fout << "###   number of iterations : " << iteration_counter << endl;
  fout << "###  max no. of iterations : " << maxiter << endl;
  if (iteration_counter >= maxiter) {
    fout << "###   max number of iterations exceeded, iteration has not "
            "converged to the desired accuracy"
         << endl;
  }
  fout << "####################################################################"
          "########################################"
          "####################################################################"
       << endl;
  fout << "### CM energy (eV) [1]   alpha_deconv (cm3/s) [2]   error_deconv "
          "(cm3/s) [3]   alpha_orig (cm3/s) [4]   "
          "error_orig (cm3/s) [5]   alpha_test (cm3/s) [6]   "
          "|1-alpha_test/alpha_orig| [7]"
       << endl;
  fout << "###-----------------------------------------------------------------"
          "----------------------------------------"
          "--------------------------------------------------------------------"
       << endl;
  for (int n = 0; n < npts; n++) {
    fout << ecm[n] << ", " << alpha_new[n] << ", " << sqrt(error_new[n]) << ", "
         << alpha_0[n] << ", " << sqrt(error_0[n]) << ", " << alpha_old[n]
         << ", " << error_old[n] << endl;
  }
  fout.close();
  cout << endl
       << endl
       << "Results written to file " << outfilename << endl
       << endl;
}
