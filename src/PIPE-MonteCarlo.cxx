/**
 * @file PIPE-MonteCarlo.cxx
 *
 * @brief Monte-Carlo simulation of molecular break-up in the PIPE experimental
setup
 *
 * For details see S. Schippers et al., ChemPhysChem 24, e202300061 (2023);
 * https://doi.org/10.1002/cphc.202300061.
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: PIPE-MonteCarlo.cxx 2039 2026-07-20 07:57:32Z iamp $
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */
#include "hydroconst.h"
#include "hydromath.h"
#include "matrix.h"
#include "buildinfo.h"
#include <chrono> // clocks and time
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional> // for "bind (random_generator, random_distribution)"
#include <iostream>
#include <random> // random number generators and distributions
#include <vector>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief KER from Coulomb explosion with intervening Auger process increasing
 * the charge of fragment 1
 *
 * @param tA time in fs when the Auger event happens (to be sampled from
 * exponential random distribution)
 * @param tch core hole lifetime in fs
 * @param r0 equilibrium distance of primary molecule in nm
 * @param q1 final charge state (after Auger process) of fragment 1 (usually the
 * heavier fragment) //
 * @param q2 charge state of fragment 2
 * @param m reduced mass in eV
 */
double KERfromCoreHoleLifetime(double tA, double tch, double r0, double q1,
                               double q2, double m) {
  const double eps = 1.0E-8;
  const int count_max = 100;
  const double esqr_four_pi_eps0 = 1.44; // in eV nm

  double KER0 = q1 * q2 * esqr_four_pi_eps0 / r0; // KER for tA=0
  if (tA < 0.01 * tch)
    return KER0;

  const double cLight = hydroconst::clight_nm_fs; // speed of light in nm/fs
  double q11 = q1 - 1.0;                   // charge before Auger emission
  double k = q11 * q2 * esqr_four_pi_eps0; // in eV nm
  if (fabs(k) < eps)
    k = eps; // to prevent later division by zero
  double t0 = r0 * sqrt(0.5 * r0 * m / k) / cLight; // in fs

  // The distance rA, where the Auger process happens, must be calculated from
  // tA. We need to sovle the equation tA = int_r0^rA dr/sqrt[2k/m(1/r0-1/r)]
  // for rA. After some algebra we arrive at tt =
  // sqrt(rr)*sqrt(rr-1.0)+ln(sqrt(rr)+sqrt(rr-1)) where tt=tA/t0  and rr =
  // rA/r0 we apply Newton's method for inverting the above equation;
  double tt = tA / t0;
  double rr = 1.0, rr_new = 1.0 + 10.0 * eps; // initial guess
  int count = 0;
  while ((count < count_max) && (fabs((rr_new - rr) / rr) > eps)) {
    count++;
    rr = rr_new;
    double srr = sqrt(rr), srr1 = sqrt(rr - 1.0);
    double f = srr * srr1 + log(srr + srr1) - tt;
    double df_drr =
        0.5 * (srr1 / srr + srr / srr1 +
               1.0 / (srr + srr1) *
                   (1.0 / srr + 1.0 / srr1)); // derivative with respect to rr
    rr_new = rr - f / df_drr;
    if (rr_new < 1.0)
      rr_new =
          1.0 + 10.0 * eps / (1.0 + count); // new inital guess (closer to 1.0)

    // cout << count << "   " <<  rr <<   "    " << rr_new  << endl;;
    // tA_new = t0*(srr*srr1+log(srr+srr1)); // check whether calcuated rA
    // produces correct tA
  }

  if (count >= count_max) {
    cout << "ERROR: Number of iterations in KERfromCoreHoleLifetime exceeded "
            "for tA = "
         << tA << " fs." << endl;
    return -1.0;
  } else {
    double rA =
        rr * r0; // internuclear distance where the Auger process happens
    double KER = q2 * esqr_four_pi_eps0 * (q11 / r0 + 1 / rA);
    // cout << tA << "  " << tA_new << "   "  << KER << endl;
    return KER;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief random polar angle
 *
 * @param r random number from a uniform distribution of real numbers from the
 * interval [0,1]
 * @param a2 ansisotropy parameter of the distribution (allowed values: -1<= a2
 * <= 2)
 *
 * @return random value for the cosine of the polar angle
 */
double random_costheta(double r, double a2, bool test_flag = false) {
  // note that in general the angle theta is calculated as a root of the cubic
  // equation r = 0.25*[2-(a2-2)*cos(theta)+a2*cos(theta)^3].
  const double eps = 1.0E-9;
  const complex<double> I(0.0, 1.0);
  double costheta = 1.0;
  if (fabs(a2) < eps) { // a2=0, isotropic distribution
    costheta = 1.0 - 2.0 * r;
    if (test_flag)
      cout << a2 << "  " << r << "  1  " << costheta << endl;
  } else if (fabs(a2 + 1.0) < eps) { // a2 = -1.0;
    complex<double> sqrtr = sqrt((1.0 - r) * r);
    double aterm = arg(1.0 - 2.0 * r + 2.0 * I * sqrtr) / 3.0;
    costheta = sqrt(3.0) * sin(aterm) - cos(aterm);
    if (test_flag)
      cout << a2 << "  " << r << "  1  " << costheta << endl;
  } else if (fabs(a2 - 2.0) < eps) { // a2 = 2.0;
    costheta = cbrt(1.0 - 2.0 * r);
    if (test_flag)
      cout << a2 << "  " << r << "   " << costheta << endl;
  } else { // general case
    double costheta1, costheta2, costheta3;
    int nsol = cuberoot(-0.25 * a2, 0.0, -0.5 + 0.25 * a2, 0.5 - r, costheta1,
                        costheta2,
                        costheta3); // returns the number of solutions
    costheta =
        (nsol == 1) ? costheta1 : costheta2; // which root to take, was looked
                                             // at by running examples
    if (test_flag)
      cout << a2 << "  " << r << "  " << nsol << "  " << costheta1 << "  "
           << costheta2 << "  " << costheta3 << endl;
  }
  if (costheta > 1.0)
    costheta = 1.0;
  if (costheta < -1.0)
    costheta = -1.0;
  return costheta;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief matrix for drift length of length L
 *
 * @param L drift length
 *
 * @return matrix with appropriately set elements
 *
 *    1, L, 0, 0, 0
 *    0, 1, 0, 0, 0
 *    0, 0 ,1, 0, 0
 *    0, 0, 0, 1, L
 *    0, 0, 0, 0, 1
 */
matrix<double> matrix_drift(double L) {
  matrix<double> m(5, 5, 0.0);
  for (unsigned int i = 0; i < 5; i++)
    m(i, i) = 1.0;
  m(0, 1) = L;
  m(3, 4) = L;
  return m;
  cout << endl << " drift: " << endl;
  for (unsigned int i = 0; i < m.get_rows(); i++) {
    for (unsigned int j = 0; j < m.get_cols(); j++)
      cout << m(i, j) << " ";
    cout << endl;
  }
  return m;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief matrix for thin lens of focal length f
 *
 * @param f focal length (make sure that f != 0.0)
 *
 * @return matrix with appropriately set elements
 *
 *    1,   0, 0,  0,   0
 *    1/f, 1, 0,  0,   0
 *    0,   0 ,1,  0,   0
 *    0,   0, 0,  1,   0
 *    0,   0, 0, -1/f, 1
 */
matrix<double> matrix_lens(double f) {
  matrix<double> m(5, 5, 0.0);
  for (unsigned int i = 0; i < 5; i++)
    m(i, i) = 1.0;
  m(1, 0) = 1.0 / f;
  m(4, 3) = -1.0 / f;
  return m;
  cout << endl << " lens: " << endl;
  for (unsigned int i = 0; i < m.get_rows(); i++) {
    for (unsigned int j = 0; j < m.get_cols(); j++)
      cout << m(i, j) << " ";
    cout << endl;
  }
  return m;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief matrix for dipole magnet
 *
 * see Wollnik, Eqs. 4.8a, 4.8b (also 4.25 for focussing)
 *
 * @param r deflection radius
 * @param phi deflection angles
 *
 * @return matrix with appropriately set elements
 *
 *      c, s*r, (1-c)*r,    0,     0
 *   -s/r,   c,       s,    0,     0
 *      0,   0,       1,    0,     0
 *      0,   0,       0,    1, r*phi
 *      0,   0,       0,    0,     1
 */
matrix<double> matrix_dipole_magnet(double r, double phi) {
  double c = cos(phi);
  double s = sin(phi);
  matrix<double> m(5, 5, 0.0);
  for (unsigned int i = 0; i < 5; i++)
    m(i, i) = 1.0;
  m(0, 0) = c;
  m(0, 1) = s * r;
  m(0, 2) = (1.0 - c) * r;
  m(1, 0) = -s / r;
  m(1, 1) = c;
  m(1, 2) = s;
  m(3, 4) = r * phi;
  return m;
  cout << endl << " magnet: " << endl;
  for (unsigned int i = 0; i < m.get_rows(); i++) {
    for (unsigned int j = 0; j < m.get_cols(); j++)
      cout << m(i, j) << " ";
    cout << endl;
  }
  return m;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Entry to Monte Carlo simulation of molecular breakup in the PIPE setup
 */
int PIPE_MonteCarlo() {

  // { for (double a2=-1.0; a2<2.05; a2+=0.1) for (double r=0.0; r <1.05;
  // r+=0.1) random_costheta(r,a2,true); } // for testing

  const double Pi = hydroconst::pi;

  const double PIPE_interact_length =
      495; // length (in mm) of interaction region
  const double PIPE_drift_demerger =
      865; // drift length (in mm) from exit of interaction region to entrance
           // of demerger
  const double PIPE_drift_detA =
      150; // drift length (in mm) from exit of interaction region to detector A
  const double PIPE_drift_detB =
      1500; // drift length (in mm) from exit of demerger to detector B
  const double PIPE_dia_interact = 6; // diameter (in mm) of the aperture at the
                                      // exit of the interaction region
  const double PIPE_dia_demerger =
      27; // diameter (in mm) of the aperture at the entrance of the demerger
  const double PIPE_dia_detectorA = 75; // diameter (in mm) of detector A
  const double PIPE_dia_detectorB = 20; // diameter (in mm) of detector B
  const double PIPE_demerger_radius =
      608.28; // bending radius (in mm) of the demerger magnet
  const double PIPE_CST_Xpos = -1500; // x position of interaction region in CST

  cout << endl
       << " Monte Carlo simulation of molecular break up in the photon-ion "
          "merged-beams setup PIPE"
       << endl;
  cout << endl
       << " The simulation yields particle distributions on two different "
          "detectors A and B.";
  cout << endl
       << " Detector A is a 2D position sensitive detector located in front of "
          "the demerging magnet.";
  cout << endl
       << " Detector B is a 1D single particle detector located behind the "
          "demerging magnet.";
  cout << endl
       << " It counts reaction products as a function of the demerging "
          "magnetic field.";
  cout << endl;
  double interact_length = PIPE_interact_length;
  double drift_demerger = PIPE_drift_demerger;
  double drift_detA = PIPE_drift_detA;
  double drift_detB = PIPE_drift_detB;
  double dia_interact = PIPE_dia_interact;
  double dia_demerger = PIPE_dia_demerger;
  double dia_detectorA = PIPE_dia_detectorA;
  double dia_detectorB = PIPE_dia_detectorB;

  char answer;
  cout << endl << " Full calculation or generation of CST input? (f/c) ....: ";
  cin >> answer;
  bool CST_input_flag = ((answer == 'c') || (answer == 'C'));

  cout << endl << " Change default geometry values? (y/n) .................: ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    cout << endl
         << " Give length of interaction region (IR) in mm ..........: ";
    cin >> interact_length;
    cout << endl
         << " Give diameter of aperture at end of IR in mm ..........: ";
    cin >> dia_interact;
    cout << endl
         << " Give drift length from IR to detector A in mm .........: ";
    cin >> drift_detA;
    cout << endl
         << " Give outer diameter of detector A in mm ...............: ";
    cin >> dia_detectorA;
    cout << endl
         << " Give drift length from IR to demerger (DEM) in mm .....: ";
    cin >> drift_demerger;
    cout << endl
         << " Give diameter of aperture at entrance of DEM ..........: ";
    cin >> dia_demerger;
    cout << endl
         << " Give drift length from DEM to detector B in mm ........: ";
    cin >> drift_detB;
    cout << endl
         << " Give outer diameter of detector B in mm ...............: ";
    cin >> dia_detectorB;
  }
  double rsqr_interact = 0.25 * dia_interact * dia_interact;
  double rsqr_detectorA = 0.25 * dia_detectorA * dia_detectorA;
  double rsqr_detectorB = 0.25 * dia_detectorB * dia_detectorB;
  double rsqr_demerger = 0.25 * dia_demerger * dia_demerger;

  // ****** input of particle properties ************
  double ion_beam_energy, frag1_mass, frag1_charge, frag2_mass,
      frag2_charge, frag_KERmin, frag_KERmax;
  double core_hole_lifetime = 1.0, equilibrium_distance = 0.1, vib_omega = 0.0;
  double frag_anisotropy; // a2 anisotropy parameter, physically real (-1..2)
  cout << endl << " Give primary ion energy in keV ........................: ";
  cin >> ion_beam_energy;
  ion_beam_energy *= 1000.0; // conversion from keV to eV
  cout << endl << " Give mass of fragment 1 in u ..........................: ";
  cin >> frag1_mass;
  cout << endl << " Give charge state of fragment 1 .......................: ";
  cin >> frag1_charge;
  cout << endl << " Give mass of fragment 2 in u ..........................: ";
  cin >> frag2_mass;
  cout << endl << " Give charge state of fragment 2 .......................: ";
  cin >> frag2_charge;
  cout << endl << " KER from core-hole lifetime? (y/n) ....................: ";
  cin >> answer;
  bool KER_from_core_hole_lifetime_flag = ((answer == 'y') || (answer == 'Y'));
  if (KER_from_core_hole_lifetime_flag) {
    cout << endl
         << " Give core-hole lifetime in fs .........................: ";
    cin >> core_hole_lifetime;
    cout << endl
         << " Give molecular equilibrium distance in nm .............: ";
    cin >> equilibrium_distance;
    cout << endl
         << " Give initial vibration frequency in eV ................: ";
    cin >> vib_omega;
  } else {
    cout << endl
         << " Give minimum KER in eV ................................: ";
    cin >> frag_KERmin;
    cout << endl
         << " Give maximum KER in eV ................................: ";
    cin >> frag_KERmax;
  }
  cout << endl << " Give anisotropy parameter a2 (-1 <= a2 <= 2) ..........: ";
  cin >> frag_anisotropy;

  double initial_mass = frag1_mass + frag2_mass;
  double initial_velocity = sqrt(2.0 * ion_beam_energy / initial_mass);

  double frag1_Bdemerger = 45.529 * sqrt(frag1_mass * ion_beam_energy * 1E-3) /
                           (PIPE_demerger_radius * 0.001 * frag1_charge) *
                           sqrt(frag1_mass / initial_mass);
  double frag2_Bdemerger = 45.529 * sqrt(frag2_mass * ion_beam_energy * 1E-3) /
                           (PIPE_demerger_radius * 0.001 * frag2_charge) *
                           sqrt(frag2_mass / initial_mass);
  double reduced_mass_eV =
      frag1_mass * frag2_mass / initial_mass * hydroconst::muc2_eV; // in eV
  double vib_width =
      vib_omega > 1E-9
          ? hydroconst::hbarc_eV_nm / sqrt(2.0 * reduced_mass_eV * vib_omega)
          : 0.0; // in nm

  // ****** input of beam properties **************

  double ion_beam_fwhm, ion_beam_divergence, ion_beam_momentum_spread,
      width_trinos_x, width_trinos_y, width_detector_slit_x;
  cout << endl << " Give FWHM of ion beam in mm ...........................: ";
  cin >> ion_beam_fwhm;
  cout << endl << " Give divergence of ion beam in mrad ...................: ";
  cin >> ion_beam_divergence;
  cout << endl << " Give momentum spread dp/p of ion beam (FWHM) ..........: ";
  cin >> ion_beam_momentum_spread;
  cout << endl << " Give width of x-Trinos slit at end of IR in mm ........: ";
  cin >> width_trinos_x;
  cout << endl << " Give width of y-Trinos slit at end of IR in mm ........: ";
  cin >> width_trinos_y;
  cout << endl << " Give width of x-slit in front of detector B in mm .....: ";
  cin >> width_detector_slit_x;

  double half_width_trinos_x = 0.5 * width_trinos_x;
  double half_width_trinos_y = 0.5 * width_trinos_y;
  double half_width_detector_slit_x = 0.5 * width_detector_slit_x;

  // ****** setup of transfer matrices ************

  double f_demerger =
      PIPE_demerger_radius /
      tan(Pi * 26.565 / 180); // focal length due to oblique entrance and exit
                              // with an angle of 26.565 deg (double focussing)
  matrix<double> m_drift_detA(matrix_drift(drift_detA));
  matrix<double> m_drift_demerger(
      matrix_drift(drift_demerger)); // from TRINOS slits to demerger
  matrix<double> m_demerger(
      matrix_lens(f_demerger) *
      matrix_dipole_magnet(PIPE_demerger_radius, 0.5 * Pi) *
      matrix_lens(f_demerger)); // demerger
  matrix<double> m_drift_detB(
      matrix_drift(drift_detB)); /// from image plane to detector B

  // *****  setup of detectors  *****
  if (CST_input_flag) {
    answer = 'r';
  } else {
    cout << endl
         << " Which detector: A, B, both A&B, or raw data? (A/B/C/r) : ";
    cin >> answer;
  }
  bool detectorB_flag = ((answer == 'b') || (answer == 'B') ||
                         (answer == 'c') || (answer == 'C'));
  bool detectorA_flag = ((answer == 'a') || (answer == 'A') ||
                         (answer == 'c') || (answer == 'C'));
  bool raw_data_flag = (!detectorB_flag) && (!detectorA_flag);

  // *****  detector A (2D detector in front of the demerging magnet)  *****
  int nbinA = 100;
  double Amin = -0.5 * dia_detectorA;
  double Adelta = dia_detectorA / (double)(nbinA);
  vector<int> detectorA1, detectorA2;
  if (detectorA_flag) {
    detectorA1.resize(nbinA * nbinA, 0); // 2D-array
    detectorA2.resize(nbinA * nbinA, 0); // 2D-array
  }
  int frag1_A_count = 0, frag2_A_count = 0;

  // ***** detector B (single particle detector behind the demerging magnet)
  // *****
  int nbinB = 1000;
  double Bdmin = -0.025,
         Bdmax = 0.025; // range of relative dB/B = dp/p variation
  double Bddelta = (Bdmax - Bdmin) / (double)(nbinB);
  vector<double> Bdemerger;
  vector<int> B1hit, B1miss, B1cut_by_demerger;
  vector<int> B2hit, B2miss, B2cut_by_demerger;
  if (detectorB_flag) {
    Bdemerger.resize(nbinB);
    B1hit.resize(nbinB, 0);
    B1miss.resize(nbinB, 0);
    B1cut_by_demerger.resize(nbinB, 0);
    B2hit.resize(nbinB, 0);
    B2miss.resize(nbinB, 0);
    B2cut_by_demerger.resize(nbinB, 0);
    for (int n = 0; n < nbinB; n++)
      Bdemerger[n] = Bdmin + (n + 0.5) * Bddelta;
  }
  int frag1_enters_demerger_count = 0, frag2_enters_demerger_count = 0;

  // ************ KER distribution ************
  vector<double> KER_eV(nbinB);
  vector<int> KER_distribution0(nbinB, 0), KER_distribution1(nbinB, 0),
      KER_distribution2(nbinB, 0);
  if (KER_from_core_hole_lifetime_flag) {
    frag_KERmax =
        KERfromCoreHoleLifetime(0.0, core_hole_lifetime, equilibrium_distance,
                                frag1_charge, frag2_charge, reduced_mass_eV);
    frag_KERmin = 0.5 * frag_KERmax * (frag1_charge - 1.0) / frag1_charge;
    frag_KERmax *= 2.0;
  }
  double frag_KERdelta = (frag_KERmax - frag_KERmin) / (double)nbinB;
  bool KER_distribution_flag = (frag_KERdelta > 0);
  for (int n = 0; KER_distribution_flag && (n < nbinB); n++) {
    KER_eV[n] = frag_KERmin + (n + 0.5) * frag_KERdelta;
  }
  //************************ initialize random generators and random
  //distributions  *******************
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937_64 randgen(
      seed); // 64bit "Mersenne Twister 19937" random number generator
  const double ln2 = log(2.0);

  uniform_real_distribution<double> dist_zero_one(0.0, 1.0);
  uniform_real_distribution<double> dist_azimut(0.0, 2.0 * Pi);
  uniform_real_distribution<double> dist_merging_section(0.0, interact_length);
  normal_distribution<double> dist_ion_beam_width(
      0.0, 0.25 * ion_beam_fwhm /
               sqrt(ln2)); // includes conversion from FWHM to Gaussian sigma
  normal_distribution<double> dist_ion_beam_momentum_spread(
      0.0, 0.25 * ion_beam_momentum_spread /
               sqrt(ln2)); // includes conversion from FWHM to Gaussian sigma
  uniform_real_distribution<double> dist_ion_beam_divergence(
      0.0, ion_beam_divergence * 0.001); // includes conversion from mrad to rad
  uniform_real_distribution<double> dist_KER(frag_KERmin, frag_KERmax);
  exponential_distribution<double> dist_core_hole_lifetime(1.0 /
                                                           core_hole_lifetime);
  normal_distribution<double> dist_mol_vibration(equilibrium_distance,
                                                 vib_width);

  int MonteCarlo_iterations;
  cout << endl << " Give the number of Monte-Carlo iterations .............: ";
  cin >> MonteCarlo_iterations;

  // ********************* open output file and write header
  // ************************
  fstream fout, fcst;
  string filename, outfilename, cstfilename, extension;

  if (CST_input_flag) {
    extension = ".pid";
  } else {
    extension = ".pmc";
  }
  cout << endl
       << " Give filename for output (" << extension
       << ") .......................: ";
  cin >> filename;

  outfilename = filename + ".pmc";
  fout.open(outfilename, fstream::out);

  if (CST_input_flag) {
    cstfilename = filename + ".pid";
    fcst.open(cstfilename, fstream::out);
  }

  // *****************************  time stamp
  // **************************************
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  fout << "####################################################################"
          "########################################"
          "###########"
       << endl;
  fout << "###  Monte Carlo simulation of molecular breakup in the PIPE setup"
       << endl;
  fout << "###" << endl;
  fout << "###            Publication : ChemPhysChem 24, e202300061 (2023)" << endl;
  fout << "###                    DOI : https://doi.org/10.1002/cphc.202300061" << endl;
  fout << "###" << endl;
  fout << "###  hydrocal revision     : " << HYDROCAL_REVISION
       << endl; // defined in buildinfo.h
  fout << "###               filename : " << outfilename << endl;
  fout << "###      start date & time : "
       << asctime(timeinfo); // endl comes with asctime
  if (detectorA_flag && detectorB_flag) {
    fout << "###              detectors : both A & B" << endl;
  } else if (detectorA_flag) {
    fout << "###              detectors : only A" << endl;
  } else if (detectorB_flag) {
    fout << "###              detectors : only B" << endl;
  } else {
    if (CST_input_flag) {
      fout << "###              detectors : as defined by CST" << endl;
    } else {
      fout << "###              detectors : none (raw data)" << endl;
    }
  }
  fout << "###  interact. length (mm) : " << interact_length << endl;
  fout << "###  apert. interact. (mm) : " << dia_interact << endl;
  fout << "###     x Trinos slit (mm) : " << width_trinos_x << endl;
  fout << "###     y Trinos slit (mm) : " << width_trinos_y << endl;
  if (!CST_input_flag) {
    fout << "###  apert. demerger  (mm) : " << dia_demerger << endl;
    fout << "### drift length detA (mm) : " << drift_detA << endl;
    fout << "###  diameter of detA (mm) : " << dia_detectorA << endl;
    fout << "###  no. of bins detect. A : " << nbinA << endl;
    fout << "### drift length dem  (mm) : " << drift_demerger << endl;
    fout << "###      demerger rho (mm) : " << PIPE_demerger_radius << endl;
    fout << "### drift length detB (mm) : " << drift_detB << endl;
    fout << "###   x detector slit (mm) : " << width_detector_slit_x << endl;
    fout << "###  diameter of detB (mm) : " << dia_detectorB << endl;
    fout << "###  no. of bins detect. B : " << nbinB << endl;
  }
  fout << "###  ion beam energy (keV) : " << ion_beam_energy * 1E-3 << endl;
  fout << "###     ion beam FWHM (mm) : " << ion_beam_fwhm << endl;
  fout << "### beam divergence (mrad) : " << ion_beam_divergence << endl;
  fout << "### momentum spread (FWHM) : " << ion_beam_momentum_spread << endl;
  fout << "###    fragment 1 mass (u) : " << frag1_mass << endl;
  fout << "###      fragment 1 charge : " << frag1_charge << endl;
  fout << "### fragm. 1 Bdemerger (G) : " << frag1_Bdemerger << endl;
  if (!CST_input_flag) {
    fout << "###    fragment 2 mass (u) : " << frag2_mass << endl;
    fout << "###      fragment 2 charge : " << frag2_charge << endl;
    fout << "### fragm. 2 Bdemerger (G) : " << frag2_Bdemerger << endl;
  }
  if (KER_from_core_hole_lifetime_flag) {
    fout << "###  core-hole lifet. (fs) : " << core_hole_lifetime << endl;
    fout << "### equilib. distance (nm) : " << equilibrium_distance << endl;
    fout << "###        hbar*omega (eV) : " << vib_omega << endl;
  } else {
    fout << "###       minimum KER (eV) : " << frag_KERmin << endl;
    fout << "###       maximum KER (eV) : " << frag_KERmax << endl;
  }
  fout << "###   anisotropy parameter : " << frag_anisotropy << endl;
  fout << "###            random seed : " << seed << endl;
  fout << "### Monte-Carlo iterations : " << MonteCarlo_iterations << endl;
  if (raw_data_flag && (!CST_input_flag)) {
    fout << "##################################################################"
            "########################################"
            "###################"
         << endl;
    fout << "###  detAxfrag1 (mm) [1]   detAyfrag1 (mm) [2]   detBxfrag1 (mm) "
            "[3]   detByfrag1 (mm) [4]   "
            "detBdeltafrag1 [5]";
    fout << "   detAxfrag2 (mm) [6]   detAyfrag2 (mm) [7]   detBxfrag2 (mm) "
            "[8]   detByfrag2 (mm) [9]   detBdeltafrag2 "
            "[10]"
         << endl;
    fout << "###---------------------------------------------------------------"
            "----------------------------------------"
            "-------------------"
         << endl;
  }

  cout << endl << endl << " now calculating, be patient ..." << endl << endl;

  double frag1_mass_SI =
      frag1_mass * hydroconst::mu_kg; // required for CST input
  double frag1_charge_SI =
      frag1_charge * hydroconst::e_As;          // required for CST input
  double sqrt_muc2 = sqrt(hydroconst::muc2_eV); // required for CST input
  int CST_counter = 0;                          // counts CST particles
  // *************************************  Monte-Carlo loop
  // ************************
  for (int i = 0; i < MonteCarlo_iterations; i++) {
    // pick a random position on the beam axis in the interaction region
    double s = dist_merging_section(randgen);

    // pick a position in the x-y-plane of the beam
    double b = dist_ion_beam_width(randgen);
    double phib = dist_azimut(randgen);
    double bx = b * cos(phib);
    double by = b * sin(phib);

    // pick angles with respect to the ion beam axis accounting for the beam
    // divergence
    double theta0 = dist_ion_beam_divergence(randgen);
    double phi0 = dist_azimut(randgen);

    // pick additional relative momentum due to the momentum spread of the ion
    // beam
    double delta0 = dist_ion_beam_momentum_spread(randgen);
    double primion_vlab =
        initial_velocity *
        (1.0 + delta0); // velocity of the primary molecular ion

    // pick a random molecular orientation
    double costheta = random_costheta(dist_zero_one(randgen), frag_anisotropy);
    double sintheta = sqrt(1.0 - costheta * costheta);
    double phi = dist_azimut(randgen);

    // pick a KER
    double frag_KER = 0;
    if (KER_from_core_hole_lifetime_flag) {
      double tA = dist_core_hole_lifetime(randgen);
      double x_vib = -1.0;
      while (x_vib < 1E-6)
        x_vib = dist_mol_vibration(randgen);
      frag_KER =
          KERfromCoreHoleLifetime(tA, core_hole_lifetime, x_vib, frag1_charge,
                                  frag2_charge, reduced_mass_eV);
    } else {
      frag_KER = dist_KER(randgen);
    }
    if (frag_KER < 0)
      continue; // negative results indicates that an error has occured (see
                // KERfromCoreHoleLifetime above)

    int nKER = 0;
    if (KER_distribution_flag) {
      nKER = floor((frag_KER - frag_KERmin) / frag_KERdelta);
      if ((nKER >= 0) && (nKER < nbinB))
        KER_distribution0[nKER]++;
    }

    // momentum conservation upon fragmentation in comoving frame: p1=-p2=|pf|
    // kinetic energy release: KER=p1*p1/(2*m1)+p2*p2/(2*m2) =
    // pf*pf/2*(m1+m2)/(m1*m2) => |pf| = sqrt[2*KER*m1*m2/(m1+m2)]
    double frag_momentum = sqrt(2.0 * frag_KER * frag1_mass * frag2_mass /
                                (frag1_mass + frag2_mass)); // in sqrt(eV*u)

    // ++++++++++ fragment 1 in the interaction region +++++++++++++++++++++++
    double frag1_vcm =
        frag_momentum /
        frag1_mass; // center-of-mass velocity of fragment 1 in sqrt(eV/u)
    double frag1_vlab_z = frag1_vcm * costheta + primion_vlab * cos(theta0);
    double frag1_vlab_x = frag1_vcm * sintheta * cos(phi) +
                          primion_vlab * sin(theta0) * cos(phi0);
    double frag1_vlab_y = frag1_vcm * sintheta * sin(phi) +
                          primion_vlab * sin(theta0) * sin(phi0);
    double frag1_time =
        s / frag1_vlab_z; // travel time to the end of the interaction region
    double x1 =
        bx +
        frag1_vlab_x *
            frag1_time; // x position at the end of the interaction region in mm
    double y1 =
        by +
        frag1_vlab_y *
            frag1_time; // y position at the end of the interaction region in mm
    double ax1 = frag1_vlab_x / frag1_vlab_z; // x angle
    double ay1 = frag1_vlab_y / frag1_vlab_z; // y angle
    double delta1 = frag1_vlab_z / initial_velocity -
                    1.0; // deviation from rigidity of test particle
    vector<double>
        v_frag1; // particle vector at the end of the interaction region
    v_frag1.push_back(x1);
    v_frag1.push_back(ax1);
    v_frag1.push_back(delta1);
    v_frag1.push_back(y1);
    v_frag1.push_back(ay1);
    bool frag1_I_flag =
        (((x1 * x1 + y1 * y1) < rsqr_interact) &&
         (fabs(x1) < half_width_trinos_x) && (fabs(y1) < half_width_trinos_y));

    //====================================================================================================================
    // CST
    if (CST_input_flag) { // write particle coordinates to CST input file, if
                          // requested
      if (frag1_I_flag) { // particle has passed the apertures at the exit of
                          // the interaction region
        CST_counter++;
        // In hydrocal, z is the ion-beam direction, x is the horizontal
        // direction, and y the vertical direction. Note that a different
        // coordinate system is used in CST. In CST, x is the ion-beam
        // direction,and y and z are the horizontal and vertical directions,
        // respectively.
        double pos_CSTx =
            PIPE_CST_Xpos *
            0.001; // in m (we are at the exit of the interaction region)
        double pos_CSTy = x1 * 0.001;                // in m
        double pos_CSTz = y1 * 0.001;                // in m
        double beta_CSTx = frag1_vlab_z / sqrt_muc2; // vx/c
        double beta_CSTy = frag1_vlab_x / sqrt_muc2; // vy/c
        double beta_CSTz = frag1_vlab_y / sqrt_muc2; // vz/c
        double mom_CSTx =
            beta_CSTx / sqrt(1.0 - beta_CSTx * beta_CSTx); // beta*gamma
        double mom_CSTy =
            beta_CSTy / sqrt(1.0 - beta_CSTy * beta_CSTy); // beta*gamma
        double mom_CSTz =
            beta_CSTz / sqrt(1.0 - beta_CSTz * beta_CSTz); // beta*gamma
        fcst << pos_CSTx << " " << pos_CSTy << " " << pos_CSTz;
        fcst << " " << mom_CSTx << " " << mom_CSTy << " " << mom_CSTz;
        fcst << " " << frag1_mass_SI << " " << frag1_charge_SI << " "
             << fabs(frag1_charge_SI) << endl;
      }
      continue; // ignore the remainder of the Monte-Carlo loop
    }
    //====================================================================================================================
    // CST

    // check whether particle hits detector A
    vector<double> v_frag1_A(m_drift_detA * v_frag1);
    bool frag1_A_flag =
        frag1_I_flag && ((v_frag1_A[0] * v_frag1_A[0] +
                          v_frag1_A[3] * v_frag1_A[3]) < rsqr_detectorA);

    // check whether particle passes entrance aperture of demerging magnet
    v_frag1 = m_drift_demerger * v_frag1; // now at entrance of demerging magnet
    bool frag1_demerger_flag =
        frag1_I_flag &&
        ((v_frag1[0] * v_frag1[0] + v_frag1[3] * v_frag1[3]) < rsqr_demerger);

    // ++++++++++ fragment 2 in the interaction region ++++++++++++++++++++++++
    double frag2_vcm =
        frag_momentum /
        frag2_mass; // center-of-mass velocity of fragment 2 in sqrt(eV/u)
    double frag2_vlab_z =
        -frag2_vcm * costheta + primion_vlab * cos(theta0);
    double frag2_vlab_x = -frag2_vcm * sintheta * cos(phi) +
                          primion_vlab * sin(theta0) * cos(phi0);
    double frag2_vlab_y = -frag2_vcm * sintheta * sin(phi) +
                          primion_vlab * sin(theta0) * sin(phi0);
    double frag2_time =
        s / frag2_vlab_z; // travel time to end of interaction region
    double x2 =
        bx + frag2_vlab_x *
                 frag2_time; // x position at the end of the interaction region
    double y2 =
        by + frag2_vlab_y *
                 frag2_time; // y position at the end of the interaction region
    double ax2 = frag2_vlab_x / frag2_vlab_z; // x angle
    double ay2 = frag2_vlab_y / frag2_vlab_z; // y angle
    double delta2 = frag2_vlab_z / initial_velocity -
                    1.0; // deviation from rigidity of test particle
    vector<double>
        v_frag2; // particle vector at the end of the interaction region
    v_frag2.push_back(x2);
    v_frag2.push_back(ax2);
    v_frag2.push_back(delta2);
    v_frag2.push_back(y2);
    v_frag2.push_back(ay2);
    bool frag2_I_flag =
        (((x2 * x2 + y2 * y2) < rsqr_interact) &&
         (fabs(x2) < half_width_trinos_x) && (fabs(y2) < half_width_trinos_y));

    // check whether particle hits detector A
    vector<double> v_frag2_A(m_drift_detA * v_frag2);
    bool frag2_A_flag =
        frag2_I_flag && ((v_frag2_A[0] * v_frag2_A[0] +
                          v_frag2_A[3] * v_frag2_A[3]) < rsqr_detectorA);

    // check whether particle passes entrance aperture of demerging magnet
    v_frag2 = m_drift_demerger * v_frag2; // now at entrance of demerging magnet
    bool frag2_demerger_flag =
        frag2_I_flag &&
        ((v_frag2[0] * v_frag2[0] + v_frag2[3] * v_frag2[3]) < rsqr_demerger);

    //++++++++++++++++++++++++++++++++++++++++++++++++++
    if (frag1_A_flag) { // fragment 1 hits detector A
      if (raw_data_flag) {
        fout << v_frag1_A[0] << ", " << v_frag1_A[3] << ", ";
      } else if (detectorA_flag) {
        int nx = floor((v_frag1_A[0] - Amin) / Adelta);
        int ny = floor((v_frag1_A[3] - Amin) / Adelta);
        detectorA1[nx * nbinA + ny]++;
      }
      frag1_A_count++;
    } else {
      if (raw_data_flag)
        fout << " --, --,";
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++
    if (frag2_A_flag) { // fragment 2 hits detector A
      if (raw_data_flag) {
        fout << v_frag2_A[0] << ", " << v_frag2_A[3] << ", ";
      } else if (detectorA_flag) {
        int nx = floor((v_frag2_A[0] - Amin) / Adelta);
        int ny = floor((v_frag2_A[3] - Amin) / Adelta);
        detectorA2[nx * nbinA + ny]++;
      }
      frag2_A_count++;
    } else {
      if (raw_data_flag)
        fout << " --, --,";
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++
    if (frag1_demerger_flag) { // fragment 1 passes the demerger entrance
                               // diaphragm
      if (raw_data_flag) {
        v_frag1 = m_drift_detB * m_demerger * v_frag1;
        fout << v_frag1[0] << ", " << v_frag1[3] << ", " << v_frag1[2] << ", ";
      } else {
        for (int n = 0; detectorB_flag && (n < nbinB);
             n++) { // this is the loop that changes the B field
          vector<double> v_frag1B(
              v_frag1); // particle coordinates at the entrance of the demerger
          v_frag1B[2] -= Bdemerger[n];      // new dB/B = dp/p
          v_frag1B = m_demerger * v_frag1B; // now at demerger exit
          if ((v_frag1B[0] * v_frag1B[0] + v_frag1B[3] * v_frag1B[3]) <
              rsqr_demerger) {
            v_frag1B = m_drift_detB * v_frag1B; // now at detector B
            if (((v_frag1B[0] * v_frag1B[0] + v_frag1B[3] * v_frag1B[3]) <
                 rsqr_detectorB) &&
                (fabs(v_frag1B[0]) < half_width_detector_slit_x)) {
              B1hit[n]++; // counts the particles that hit the detector
              if (KER_distribution_flag && (nKER >= 0) && (nKER < nbinB))
                KER_distribution1[nKER]++;
            } else {
              B1miss[n]++; // counts the particles that do not hit the detector
            }
          } else {
            B1cut_by_demerger[n]++; // counts the particles that do not pass the
                                    // demerger exit
          }
        } // end for(n...)
      }
      frag1_enters_demerger_count++;
    } else {
      if (raw_data_flag)
        fout << " -- ,--, --,";
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++
    if (frag2_demerger_flag) { // fragment 2 passes the demerger entrance
                               // diaphragm
      if (raw_data_flag) {
        v_frag2 = m_drift_detB * m_demerger * v_frag2;
        fout << v_frag2[0] << ", " << v_frag2[3] << ", " << v_frag2[2];
      } else {
        for (int n = 0; detectorB_flag && (n < nbinB);
             n++) { // this is the loop that changes the B field
          vector<double> v_frag2B(
              v_frag2); // particle coordinates at the entrance of the demerger
          v_frag2B[2] -= Bdemerger[n];      // new dB/B = dp/p
          v_frag2B = m_demerger * v_frag2B; // now at demerger exit
          if ((v_frag2B[0] * v_frag2B[0] + v_frag2B[3] * v_frag2B[3]) <
              rsqr_demerger) {
            v_frag2B = m_drift_detB * v_frag2B; // now at detector B
            if (((v_frag2B[0] * v_frag2B[0] + v_frag2B[3] * v_frag2B[3]) <
                 rsqr_detectorB) &&
                (fabs(v_frag2B[0]) < half_width_detector_slit_x)) {
              B2hit[n]++; // counts the particles that hit the detector
              if (KER_distribution_flag && (nKER >= 0) && (nKER < nbinB))
                KER_distribution2[nKER]++;
            } else {
              B2miss[n]++; // counts the particles that do not hit the detector
            }
          } else {
            B2cut_by_demerger[n]++; // counts the particles that do not pass the
                                    // demerger exit
          }
        } // end for(n...)
      }
      frag2_enters_demerger_count++;
    } else {
      if (raw_data_flag)
        fout << " -- ,--, --";
    }
    if (raw_data_flag)
      fout << endl;
  } // end for(i ...) end of Monte-Carlo loop
  ///////////////////////////////////////////////////////////////////////////////////////
  if (CST_input_flag) {
    fcst.close();
    fout << "###   no. of CST particles : " << CST_counter << endl;
    fout << "###         CST input file : " << cstfilename << endl;
    fout << "##################################################################"
            "########################################"
            "#############"
         << endl;
    fout.close();
    cout << " " << CST_counter
         << " CST particle coordinates written to CST input file "
         << cstfilename << "." << endl;
    cout << endl
         << " Monte-Carlo input parameters written to file: " << outfilename
         << endl
         << endl;
    return 0;
  }

  fout << "####################################################################"
          "########################################"
          "###########"
       << endl;
  fout << "###  fragm. 1 on detect. A : " << frag1_A_count << endl;
  fout << "###  fragm. 2 on detect. A : " << frag2_A_count << endl;
  fout << "###  fragm. 1 into demerger: " << frag1_enters_demerger_count
       << endl;
  fout << "###  fragm. 2 into demerger: " << frag2_enters_demerger_count
       << endl;

  if (raw_data_flag) {
    fout.close();
    cout << endl
         << endl
         << " Monte-Carlo raw data written to output file: " << outfilename
         << endl
         << endl;
    return 0;
  }

  fstream fdetector;
  //////////////////////////////////////////////////////////////////////////////////////
  // output of fragment 1 on detector A
  if (detectorA_flag) {
    outfilename = filename + ".pmc1";
    fdetector.open(outfilename, fstream::out);
    for (int i = 0; i < nbinA; i++) {
      for (int j = 0; j < nbinA - 1; j++) {
        fdetector << detectorA1[i * nbinA + j] << ", ";
      }
      fdetector << detectorA1[i * nbinA + nbinA - 1] << endl;
    }
    fdetector.close();
    cout << " Data for fragment 1 on detector A stored in " << outfilename
         << endl;
  }
  //////////////////////////////////////////////////////////////////////////////////////
  // output of fragment 2 on detector A
  if (detectorA_flag) {
    outfilename = filename + ".pmc2";
    fdetector.open(outfilename, fstream::out);
    for (int i = 0; i < nbinA; i++) {
      for (int j = 0; j < nbinA - 1; j++) {
        fdetector << detectorA2[i * nbinA + j] << ", ";
      }
      fdetector << detectorA2[i * nbinA + nbinA - 1] << endl;
    }
    fdetector.close();
    cout << " Data for fragment 2 on detector A stored in " << outfilename
         << endl;
  }
  //////////////////////////////////////////////////////////////////////////////////////
  // detector B output
  if (detectorB_flag) {
    fout << "##################################################################"
            "########################################"
            "##################################################################"
            "################"
         << endl;
    fout << "###  Delta_B1/B1 (G) [1]    hitB1 [2]   missB1 [3]   cut1dem [4]  "
            " Delta_B2/B2 (G) [5]   hitB2 [6]   "
            "missB2 [7]   cut2dem [8]   KER (eV) [9]   KERdistr0 [10]   "
            "KERdistr1 [11]   KERdistr2 [12]"
         << endl;
    fout << "###---------------------------------------------------------------"
            "----------------------------------------"
            "------------------------------------------------------------------"
            "-------------------"
         << endl;
    for (int n = 0; n < nbinB; n++) {
      fout << Bdemerger[n] << ", " << B1hit[n] << ", " << B1miss[n] << ", "
           << B1cut_by_demerger[n] << ", ";
      fout << Bdemerger[n] << ", " << B2hit[n] << ", " << B2miss[n] << ", "
           << B2cut_by_demerger[n];
      if (KER_distribution_flag) {
        fout << ", " << KER_eV[n] << ", " << nbinB * KER_distribution0[n]
             << ", " << KER_distribution1[n] << ", " << KER_distribution2[n]
             << endl;
      } else {
        fout << ", --, --, --, --" << endl;
      }
    }
    cout << " Summary and detector-B data written to output file " << filename
         << ".pmc." << endl
         << endl;
  } else {
    cout << " Summary written to output file " << filename << ".pmc." << endl
         << endl;
  }
  fout.close();

  return 0;
}
