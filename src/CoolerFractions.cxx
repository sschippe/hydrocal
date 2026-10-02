/**
 * @file CoolerFractions.cxx
 *
 * @brief Surviving fractions of nl-Rydberg states populated in storage-ring
recombination experiments
 *
 * For details see Schippers et al., Astrophys. J 555, 1027 (2001);
 * https://doi.org/10.1086/321512.
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: CoolerFractions.cxx 2039 2026-07-20 07:57:32Z iamp $
// SPDX-License-Identifier: MIT
 @endverbatim
 *
 */

#include "RRDRxsec.h"
#include "fieldion.h"
#include "hydroconst.h"
#include "hydromath.h"
#include "matrix.h"
#include "radrate.h"
#include "buildinfo.h"
#include <chrono> // clocks and time
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using hydroconst::pi;
static const double clight = hydroconst::clight_cm_s;

// #define printcasc 1

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Classical cut-off quantum number for a hydrogenic ion in an electric
 * field
 *
 * @param efield electric field strength in V/cm
 * @param z nuclear charge
 */
int nf(double efield, double z) {
  return int(sqrt(sqrt(Fau * z * z * z / efield / 9)) + 0.5);
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Mean cut-off quntum number for a hydrogenic ion in an electric field
 *
 * Formula developed by Andreas Wolf, Habilitation Thesis, University of
 * Heidelberg, 1992, page 75.
 *
 * @param tau ion-flight time (in s) from cooler to field-ionizing field
 * @param z nuclear charge
 * @param nf calssical cut-off quantum number (used as start value in iteration)
 */
int ngamma(double tau, double z, int nf) {
  double fac = pow(z, 4) * tau * 0.357 * 2.142e10;
  double n_new = nf, n_old = 0;
  while (fabs(n_new - n_old) > 1.E-4) {
    n_old = n_new;
    n_new = fac * exp(-1.36 * log(n_old - 1)) * (0.481 + 0.68 * log(n_old - 1));
    n_new = exp(log(n_new) / 3);
  }
  return int(n_new);
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief writes input data in human readable form to a *.txt file
 */
void write_inputdata_to_file(ofstream &ftxt, const string fnroot, const string storage_ring,
                             double qion, double mion, double ecool,
                             double bpar, int nmax, int ncasc, int cutoff,
                             int ndip = 0) {
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  ftxt << "#################################################################################################" << endl;
  ftxt << "### Input data for calculation of subshell dependent field ionization" << endl;
  ftxt << "###" << endl;
  ftxt << "###            Publication : Astrophys. J 555, 1027 (2001)" << endl;
  ftxt << "###                    DOI : https://doi.org/10.1086/321512" << endl;
  ftxt << "###" << endl;
  ftxt << "###  hydrocal revision     : " << HYDROCAL_REVISION << endl; // defined in buildinfo.h
  ftxt << "###               filename : " << fnroot << endl;
  ftxt << "###      start date & time : " << asctime(timeinfo);
  ftxt << "###           storage ring : " << storage_ring << endl;
  if (ndip > 0) {
    ftxt << "###      number of dipoles : " << ndip << endl;
  }
  ftxt << "###             ion charge : " << qion << endl;
  ftxt << "###           ion mass (u) : " << mion << endl;
  ftxt << "###    cooling energy (eV) : " << ecool << endl;
  ftxt << "###    guiding field  (mT) : " << bpar << endl;
  ftxt << "###                  max n : " << nmax << endl;
  if (ncasc==-1) {
    ftxt << "###          cascade depth : full via matrix exponentiation" << endl;
  }
  else if (ncasc<-1) {
    ftxt << "###          cascade depth : full via recursion" << endl;
  }
  else {
    ftxt << "###          cascade depth : " << ncasc << endl;
  }
  if (cutoff>0) {
    ftxt << "###                 cutoff : hard" << endl;
  } else {
    ftxt << "###                 cutoff : soft" << endl;
  }
  ftxt << "#################################################################################################" << endl;
  ftxt << "###          batch command : hydrocal < " << fnroot << ".hcin &" << endl;
  ftxt << "###" << endl;
  ftxt << "### The following output files (among others) will be created: " << endl;
  ftxt << "###   " << fnroot<< ".fn  : n-resolved survival probabilities averaged over l" << endl;
  ftxt << "###   " << fnroot<< ".fnl : nl-resolved survival probabilities for input to autostructure"  << endl;
  ftxt << "###   " << fnroot<< ".fm1 : nl-resolved survival probabilities in matrix form for 2D plotting" << endl; 
  ftxt << "#################################################################################################" << endl;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief writes full specification of field ionizing magnet to hydrocal input
 * file for field-ionization option 3)
 *
 * @param qion charge state of primary ion
 * @param mion_u mass of primary ion in u
 * @param ecool_eV cooling energy in eV
 * @param z1_cm distance from beginning of previuos device to magnet entrance in
 * cm
 * @param z2_cm distance from end of previous device to magnet entrance in cm
 * (differs from z1 only if previous device is the cooler)
 * @param dz_cm length of magnet in ion-beam direction in cm
 * @param r_cm bending radius in cm
 * @param phi_deg deflection angle in deg (a negative value signals an
 * electrostatic deflector)
 * @param Ef_V_cm electric field strength (will be calculated internally if set
 * to zero)
 * @param nF field ionization cut-off quantum number (nf<0 signals soft cutoff)
 * @param nmax maximum n to be considered (usually higher than highest nF)
 * @param ncasc number of cascade steps (usually set to 0)
 * @param fout output file (to be used as hydrocal input file in the next
 * @param ftxt text output file for file header 
 * calculation step)
 */
void field_ionizing_magnet(double qion, double mion_u, double ecool_eV,
                           double z1_cm, double z2_cm, double dz_cm,
                           double r_cm, double phi_deg, double Ef_V_cm, int nF,
                           int nmax, int ncasc, const string &fn,
                           const string &fn_old, ostream &fout, ofstream &ftxt) {
  fout << "3\n";
  fout << setw(2) << qion << "\n";
  fout << setw(6) << fixed << setprecision(2) << mion_u << "\n";
  fout << setw(9) << fixed << setprecision(2) << ecool_eV << "\n";
  fout << setw(5) << fixed << setprecision(1) << z1_cm << " " << setw(5)
       << fixed << setprecision(1) << z2_cm << " " << setw(6) << fixed
       << setprecision(2) << r_cm << "\n";
  fout << setw(6) << fixed << setprecision(2) << phi_deg << "\n";
  if (Ef_V_cm > 1.0E-6) { // use Efield as specified
    fout << "y\n";
    fout << setw(5) << fixed << setprecision(1) << dz_cm << "\n";
    fout << "y\n";
    fout << setw(9) << fixed << setprecision(2) << Ef_V_cm << "\n";
  } else { // compute Efield internally
    fout << "n\n";
    fout << "n\n";
  }
  ftxt << "### electric field (kV/cm) : " << fabs(Ef_V_cm)*1.0E-3 << endl;
  ftxt << "### field ioniz. cut off n : " << abs(nF) <<endl;
  fout << setw(4) << nF << "\n";
  fout << setw(4) << nmax << "\n";
  fout << "y\n";
  fout << setw(4) << ncasc << "\n";
  if (fn_old.empty()) {
    fout << "n\n";
  } else {
    fout << "y\n";
    fout << fn_old << "\n";
  }
  fout << fn << "\n";
  fout << "n\n";

  ftxt <<"##################################################" << endl;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at CSR cooler
 *
 * Interactive generation of a hydrocal batch-file for the calculation of
 * field-ionization survival probabilities behind the CSR electron cooler
 */
void setup_batch_csr(void) {
  const double melectron = hydroconst::mec2_eV;
  const double atomic_mass_unit = hydroconst::muc2_eV;

  // geometric quantities, all lengths are in cm
  // some values taken from the PhD-Thesis of Sunny Saurabh (Univ. Heidelberg
  // 2019, doi: 10.11588/heidok.00026791)
  const double dz_c = 110.0; // length of the cooler
  const double z_c_dc =
      55.3; // distance from the center of the cooler to the demerging coils
  const double dz_dc = 29.0;  // length of demerging coils in ion-beam direction
  const double Bfac_dc = 0.5; // B-field ratio demerger/cooler

  const double z_c_cc1 = 85.0; // distance from the center of the cooler to the
                               // 1st compensator coils
  const double dz_cc1 =
      11.5; // length of 1st compensator coils in ion-beam direction
  const double Bfac_cc1 = 2.2; // B-field ratio compensator1/cooler

  const double z_c_cc2 = 105.0; // distance from the center of the cooler to the
                                // 2nd compensator coils
  const double dz_cc2 =
      11.5; // length of 2nd compensator coils in ion-beam direction
  const double Bfac_cc2 = 1.5; // B-field ratio compensator2/cooler

  const double z_c_06 =
      215.5; // distance from the center of the cooler to the 6-deg deflector
  const double dz_06 = 20.3;      // length of orbit in 06-deg deflector
  const double rho_06 = 200.0;    // bending radius of 6-deg deflector;
  const double EperkeV_06 = 10.0; // electrical field in V/cm per keV ion energy
                                  // (36 kV/12cm at 300 kV acceleration voltage)

  const double z_06_39 = 110.7;   // distance between 6-deg and 39-deg deflector
  const double dz_39 = 67.7;      // length of orbit in 39-deg deflector
  const double rho_39 = 100.0;    // bending radius of 39-deg deflector;
  const double EperkeV_39 = 20.0; // electrical field in V/cm per keV ion energy
                                  // (36 kV/6cm at 300 kV acceleration voltage)

  bool defl39_flag = false; // set true if field ionization data from the 39-deg
                            // deflector are required

  double qion, mion, eion, Bcooler;
  printf("\n");
  printf("\n Setup of a batch file for the calculation of nl-specific "
         "detection probabilities at the CSR electron cooler");
  printf("\n");
  printf("\n Give ion charge .........................: ");
  scanf("%lf", &qion);
  printf("\n Give ion mass in u ......................: ");
  scanf("%lf", &mion);
  printf("\n Give ion energy in keV ..................: ");
  scanf("%lf", &eion);
  eion *= 1000.0; // conversion of keV to eV
  printf("\n Give magnetic field in the cooler in mT .: ");
  scanf("%lf", &Bcooler);
  Bcooler *= 1E-7; // conversion of mT to Vs/(cm)^2

  double ecool = eion * melectron / (mion * atomic_mass_unit);
  double gamma = 1.0 + ecool / melectron;
  double vion = clight * sqrt(1.0 - 1.0 / (gamma * gamma)); // in cm/s

  printf("\n The relativistic gamma factor is: %9.8f\n", gamma);
  printf("\n The ion velocity is ............: %9.5g cm/s\n", vion);
  printf("\n The cooling energy is ..........: %9.5g eV\n", ecool);

  printf("\n The flight times are (in the ion's frame of reference)");

  printf("\n through the cooler ............................: %6.2f mus",
         dz_c * 1e6 / vion / gamma);
  printf("\n through the demerging coils ...................: %6.2f mus",
         dz_dc * 1e6 / vion / gamma);
  printf("\n through the 1st compensator coils .............: %6.2f mus",
         dz_cc1 * 1e6 / vion / gamma);
  printf("\n through the 2nd compensator coils .............: %6.2f mus",
         dz_cc2 * 1e6 / vion / gamma);
  printf("\n through the  6-deg deflector ..................: %6.2f mus",
         dz_06 * 1e6 / vion / gamma);
  if (defl39_flag) {
    printf("\n through the 39-deg deflector ..................: %6.2f mus",
           dz_39 * 1e6 / vion / gamma);
  }
  double dt_dc = z_c_dc / vion / gamma;
  printf("\n from the cooler center to the demerging coils .: %6.2f mus",
         dt_dc * 1e6);
  double dt_cc1 = z_c_cc1 / vion / gamma;
  printf("\n from the cooler center to the 1st comp. coils .: %6.2f mus",
         dt_cc1 * 1e6);
  double dt_cc2 = z_c_cc2 / vion / gamma;
  printf("\n from the cooler center to the 2nd comp. coils .: %6.2f mus",
         dt_cc2 * 1e6);
  double dt_06 = z_c_06 / vion / gamma;
  printf("\n from the cooler center to the  6-deg deflector.: %6.2f mus",
         dt_06 * 1e6);
  if (defl39_flag) {
    double dt_39 = (z_c_06 + dz_06 + z_06_39) / vion / gamma;
    printf("\n from the cooler center to the 39-deg deflector.: %6.2f mus",
           dt_39 * 1e6);
  }
  // bending radii from p*c = zeff*c*B*rho
  double rho =
      sqrt(2.0 * mion * eion * atomic_mass_unit) / (qion * clight * Bcooler);
  double rho_dc = rho / Bfac_dc;
  double rho_cc1 = rho / Bfac_cc1;
  double rho_cc2 = rho / Bfac_cc2;

  // deflection angles
  double phi_dc = asin(dz_dc / rho_dc) * 180.0 / pi;
  double phi_cc2 = asin(dz_cc2 / rho_cc2) * 180.0 / pi;
  double phi_cc1 = phi_dc + phi_cc2;
  printf("\n\n The bending radii and deflection angles are ");
  printf("\n of the demerging coils ..................: %7.1f cm, %6.3f deg",
         rho_dc, phi_dc);
  printf("\n of the 1st compensating coils ...........: %7.1f cm, %6.3f deg",
         rho_cc1, phi_cc1);
  printf("\n of the 2nd compensating coils ...........: %7.1f cm, %6.3f deg",
         rho_cc2, phi_cc2);
  printf("\n of the 6-deg deflector...................: %7.1f cm, %6.3f deg",
         rho_06, 6.0);
  if (defl39_flag) {
    printf("\n of the 39-deg deflector..................: %7.2f cm, %6.3f deg",
           rho_39, 39.0);
  }
  // motional electric fields
  double ef_dc = Bcooler * Bfac_dc * vion;
  double ef_cc1 = Bcooler * Bfac_cc1 * vion;
  double ef_cc2 = Bcooler * Bfac_cc2 * vion;
  double ef_06 = EperkeV_06 * eion * 0.001 / qion;
  double ef_39 = EperkeV_39 * eion * 0.001 / qion;
  // cut-off quantum numbers (lowest n ionized)
  int nF_dc = nf(ef_dc, qion);
  int nF_cc1 = nf(ef_cc1, qion);
  int nF_cc2 = nf(ef_cc2, qion);
  int nF_06 = nf(ef_06, qion);
  int nF_39 = nf(ef_39, qion);
  printf("\n\n The motional electric fields and cut-off quantum numbers are ");
  printf("\n in the demerging coils ..................: %9.2f kV/cm",
         ef_dc * 1E-3);
  printf(", nF = %3d", nF_dc);
  printf("\n in the 1st compensator coils ............: %9.2f kV/cm",
         ef_cc1 * 1E-3);
  printf(", nF = %3d", nF_cc1);
  printf("\n in the 2nd compensator coils ............: %9.2f kV/cm",
         ef_cc2 * 1E-3);
  printf(", nF = %3d", nF_cc2);
  printf("\n in the  6-deg deflector .................: %9.2f kV/cm",
         ef_06 * 1E-3);
  printf(", nF = %3d", nF_06);
  if (defl39_flag) {
    printf("\n in the 39-deg deflector .................: %9.2f kV/cm",
           ef_39 * 1E-3);
    printf(", nF = %3d", nF_39);
  }
  int nmax, ncasc;
  string fn, fn_old, fnroot;
  char answer;
  cout << "\n\n Give maximum main quantum number ........: ";
  cin >> nmax;
  cout << "\n Give number of cascade steps ............: ";
  cin >> ncasc;
  if (ncasc>nmax) ncasc=nmax;
  cout << "\n Hard or soft cut-off ? (h/s) ............: ";
  cin >> answer;
  if ((answer == 's') || (answer == 'S')) {
     // nF<=0 signals soft cut-off
    nF_dc  = -abs(nF_dc);
    nF_cc1 = -abs(nF_cc1);
    nF_cc2 = -abs(nF_cc2);
    nF_06  = -abs(nF_06);
    nF_39  = -abs(nF_39);
  }
  cout << "\n Give filename for output (*.fnl) ........: ";
  cin >> fnroot;

  fn = fnroot + ".hcin";

  ofstream fout(fn);
  ofstream ftxt(fnroot+".txt");
  write_inputdata_to_file(ftxt, fnroot, "CSR", qion, mion, ecool, Bcooler * 1E7, nmax,
                          ncasc, nF_39);
  
  fout << "4\n";
  fn_old.clear();
  ftxt << "### cooler toriod" << endl;
  fn = fnroot + "_dc";
  field_ionizing_magnet(qion, mion, ecool, z_c_dc + 0.5 * dz_c,
                        z_c_dc - 0.5 * dz_c, dz_dc, rho_dc, phi_dc, ef_dc,
                        nF_dc, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### 1st compensator" << endl;
  double z12 = z_c_cc1 - z_c_dc; // distance from demerger to 1st compensator
  fn_old = fn;
  fn = fnroot + "_c1";
  field_ionizing_magnet(qion, mion, ecool, z12, z12, dz_cc1, rho_cc1, phi_cc1,
                        ef_cc1, nF_cc1, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### 2nd compensator" << endl;
  z12 = z_c_cc2 - z_c_cc1; // distance from 1st to 2nd compensator
  fn_old = fn;
  fn = fnroot + "_c2";
  field_ionizing_magnet(qion, mion, ecool, z12, z12, dz_cc2, rho_cc2, phi_cc2,
                        ef_cc2, nF_cc2, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### 6-deg deflector" << endl;
  z12 = z_c_06 - z_c_cc2; // distance from 2nd compensator to 6-deg deflector
  fn_old = fn;
  fn = fnroot;
  field_ionizing_magnet(qion, mion, ecool, z12, z12, dz_06, rho_06, -6.0, ef_06,
                        nF_06, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### 39-deg deflector" << endl;
  if (defl39_flag) {
    z12 = z_06_39 + dz_06; // distance from 6-deg to 39-deg deflector
    fn_old = fn;
    fn = fnroot + "_39";
    field_ionizing_magnet(qion, mion, ecool, z12, z12, dz_39, rho_39, -39.0,
                          ef_39, nF_39, nmax, ncasc, fn, fn_old, fout, ftxt);
  }

  fout << "0\n";
  fout << "0\n";

  fout.close();
  ftxt.close();

  cout << "\n\n run batch job by issuing the command: nice hydrocal <" << fnroot
       << ".hcin >" << fnroot << ".log &\n";
  cout << "--------------------------------------------------------------------"
          "------------------------------\n\n";
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at CSRm cooler
 *
 * Interactive generation of a hydrocal batch-file for the calculation of
 * field-ionization survival probabilities behind the CSRm electron cooler
 */
void setup_batch_csrm(void) {
  const double melectron = hydroconst::mec2_eV;

  // geometric quantities, all lengths are in cm
  const double z1_t = 400.0;  // distance from beginning of cooler to toroid
  const double z2_t = 10.0;   // distance from end of cooler to toroid
  const double dz_t = 90.0;   // distance within toroid
  const double z_d = 800.0;   // distance from toroid to dipole
  const double rho_d = 760.0; // bending radius of dipole magnet;
  const double phi_d = 22.5;  // bending angle of dipole magnet;

  double qion, mion, ecool, bpar;
  printf("\n");
  printf("\n Setup of a batch file for the calculation of nl-specific "
         "detection probabilities at the CSRm electron cooler");
  printf("\n");
  printf("\n Give ion charge .........................: ");
  scanf("%lf", &qion);
  printf("\n Give ion mass in u ......................: ");
  scanf("%lf", &mion);
  printf("\n Give cooling energy in eV ...............: ");
  scanf("%lf", &ecool);
  printf("\n Give magnetic guiding field in mT .......: ");
  scanf("%lf", &bpar);

  bpar *= 1E-7; // conversion of mT to Vs/(cm)^2
  double gamma = 1.0 + ecool / melectron;
  double vion = clight * sqrt(1.0 - 1.0 / (gamma * gamma)); // in cm/s
  double ef_t = 0.5 * bpar * vion; // E field in toroid n V/cm

  double ef_d =
      1.036427E-12 * gamma * vion * vion * mion / qion / rho_d; // in Dipole

  // calculate the length of the recombined ion's path through the dipole magnet
  double dz_d = rho_d * qion / (qion - 1.0) *
                asin((qion - 1.0) / qion * cos(phi_d * pi / 180.0));

  printf("\n The relativistic gamma factor is: %9.8f\n", gamma);
  printf("\n The ion velocity is ............: %9.5g cm/s\n", vion);

  printf("\n The flight times are (in the ion's frame of reference)");

  printf("\n from the beginning of the cooler to the toroid : %6.2f ns",
         z1_t * 1e9 / vion / gamma);
  printf("\n from the end of the cooler to the toroid ......: %6.2f ns",
         z2_t * 1e9 / vion / gamma);
  printf("\n from the toroid to the dipole .................: %6.2f ns",
         z_d * 1e9 / vion / gamma);
  printf("\n through the cooler ............................: %6.2f ns",
         (z1_t - z2_t) * 1e9 / vion / gamma);
  printf("\n through the toroid ............................: %6.2f ns",
         dz_t * 1e9 / vion / gamma);
  printf("\n through the dipole ............................: %6.2f ns",
         dz_d * 1e9 / vion / gamma);

  double dt_t, dt_d;
  double zsum = 0.5 * (z1_t + z2_t);
  dt_t = zsum / vion / gamma;
  printf("\n from the cooler center to the toroid ..........: %6.2f ns",
         dt_t * 1e9);
  zsum += z_d;
  dt_d = zsum / vion / gamma;
  printf("\n from the cooler center to the dipole ..........: %6.2f ns",
         dt_d * 1e9);

  int nF_t = nf(ef_t, qion);
  int nF_d = nf(ef_d, qion);
  printf("\n\n The motional electric fields and cut-off quantum numbers are ");
  printf("\n in the toroid ..................: %9.2f kV/cm", ef_t * 1E-3);
  printf(", nF = %3d", nF_t);
  printf("\n in the charge analyzing dipole .: %9.2f kV/cm", ef_d * 1E-3);
  printf(", nF = %3d", nF_d);

  int nmax, ncasc;
  char answer;
  printf("\n\n Give maximum main quantum number ........: ");
  scanf("%d", &nmax);
  printf("\n Give number of cascade steps ............: ");
  scanf("%d", &ncasc);
  printf("\n Hard or soft cut-off ? (h/s) ............: ");
  scanf(" %c", &answer);
  if ((answer == 's') || (answer == 'S')) {
    // nF <= 0 signals soft cut-off
    nF_d = -abs(nF_d);
    nF_t = -abs(nF_t);
  }

  string fn, fn_old, fnroot;
  cout << "\n Give filename for output (*.fnl) ........: ";
  cin >> fnroot;

  fn = fnroot + ".hcin";
  ofstream fout(fn);
  ofstream ftxt(fnroot+".txt");
  write_inputdata_to_file(ftxt, fnroot, "CSRm", qion, mion, ecool, bpar * 1E7, nmax,
                          ncasc, nF_t);

  fout << "4\n";

  ftxt << "### cooler toriod" << endl;
  fn_old.clear();
  fn = fnroot + "_t";
  field_ionizing_magnet(qion, mion, ecool, z1_t, z2_t, dz_t, 100.0, 45.0, ef_t,
                        nF_t, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### analyzing dipole magnet" << endl;
  fn_old = fn;
  fn = fnroot + "_d";
  field_ionizing_magnet(qion, mion, ecool, z_d, z_d, dz_d, rho_d, phi_d, -ef_d,
                        nF_d, nmax, ncasc, fn, fn_old, fout, ftxt);


  fout << "0\n";
  fout << "0\n";
  fout.close();
  ftxt.close();

  cout << "\n\n run batch job by issuing the command: nice hydrocal <" << fnroot
       << ".hcin >" << fnroot << ".log &\n";
  cout << "--------------------------------------------------------------------"
          "------------------------------\n\n";
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at TSR cooler
 *
 * Interactive generation of a hydrocal batch-file for the calculation ofx1
 * field-ionization survival probabilities behind the TSR electron cooler
 */
void setup_batch_tsrc(void) {
  const double melectron = hydroconst::mec2_eV;

  const double rho = 115.0;  // bending radius of dipole magnet;
  const double phi = 45.0;   // bending angle of dipole magnet;
  const double z1_t = 170.0; // distance from beginning of cooler to toroid
  const double z2_t = 20.0;  // distance from end of cooler to toroid
  const double dz_t = 50.0;  // 56.5; // distance within toroid
  const double z_c1 = 82.5;  // distance from toroid to KDX1
  const double dz_c1 = 34.0; // distance within KDX1
  const double z_c2 = 50.5;  // distance from KDX1 to KDX2
  const double dz_c2 = 18.0; // distance within KDX2
  const double z_d = 244.0;  // distance from KDX2 to dipole

  double qion, mion, ecool, bpar;

  printf("\n");
  printf("\n Setup of a batch file for the calculation of nl-specific "
         "detection probabilities at the TSR electron cooler");
  printf("\n");
  printf("\n Give ion charge .........................: ");
  scanf("%lf", &qion);
  printf("\n Give ion mass in u ......................: ");
  scanf("%lf", &mion);
  printf("\n Give cooling energy in eV ...............: ");
  scanf("%lf", &ecool);
  printf("\n Give magnetic guiding field in mT .......: ");
  scanf("%lf", &bpar);

  bpar *= 1E-7; // conversion of mT to Vs/(cm)^2
  double gamma = 1.0 + ecool / melectron;
  double vion = clight * sqrt(1.0 - 1.0 / (gamma * gamma)); // in cm/s
  double efield = bpar * vion;                              // in V/cm

  double ef_t = 0.5 * efield; // 0.36*efield
  double ef_c1 = 0.8 * efield;
  double ef_c2 = 1.6 * efield;
  efield = 1.036427E-12 * gamma * vion * vion * mion / qion / rho;

  // calculate the length of the recombined ion's path through the dipole magnet
  double dz_d = rho * qion / (qion - 1.0) *
                asin((qion - 1.0) / qion * cos(phi * pi / 180.0));

  printf("\n The relativistic gamma factor is: %9.8f\n", gamma);
  printf("\n The ion velocity is ............: %9.5g cm/s\n", vion);

  printf("\n The flight times are (in the ion's frame of reference)");

  printf("\n from the beginning of the cooler to the toroid : %6.2f ns",
         z1_t * 1e9 / vion / gamma);
  printf("\n from the end of the cooler to the toroid ......: %6.2f ns",
         z2_t * 1e9 / vion / gamma);
  printf("\n from the toroid to the first correction magnet.: %6.2f ns",
         z_c1 * 1e9 / vion / gamma);
  printf("\n from the first to the second correction magnet.: %6.2f ns",
         z_c2 * 1e9 / vion / gamma);
  printf("\n from the second correction magnet to the dipole: %6.2f ns",
         z_d * 1e9 / vion / gamma);
  printf("\n through the cooler ............................: %6.2f ns",
         (z1_t - z2_t) * 1e9 / vion / gamma);
  printf("\n through the toroid ............................: %6.2f ns",
         dz_t * 1e9 / vion / gamma);
  printf("\n through the first correction magnet............: %6.2f ns",
         dz_c1 * 1e9 / vion / gamma);
  printf("\n through the second correction magnet...........: %6.2f ns",
         dz_c2 * 1e9 / vion / gamma);
  printf("\n through the dipole ............................: %6.2f ns",
         dz_d * 1e9 / vion / gamma);

  double dt_t, dt_c1, dt_c2, dt_d;
  double zsum = 0.5 * (z1_t + z2_t);
  dt_t = zsum / vion / gamma;
  printf("\n from the cooler center to the toroid ..........: %6.2f ns",
         dt_t * 1e9);
  zsum += z_c1;
  dt_c1 = zsum / vion / gamma;
  printf("\n from the cooler center to the 1st correction m.: %6.2f ns",
         dt_c1 * 1e9);
  zsum += z_c2;
  dt_c2 = zsum / vion / gamma;
  printf("\n from the cooler center to the 2nd correction m.: %6.2f ns",
         dt_c2 * 1e9);
  zsum += z_d;
  dt_d = zsum / vion / gamma;
  printf("\n from the cooler center to the dipole ..........: %6.2f ns",
         dt_d * 1e9);

  int nF_t = nf(ef_t, qion);
  int nF_c1 = nf(ef_c1, qion);
  int nF_c2 = nf(ef_c2, qion);
  int nF_d = nf(efield, qion);
  printf("\n\n The motional electric fields and cut-off quantum numbers are ");
  printf("\n in the toroid ..................: %9.2f kV/cm", ef_t * 1E-3);
  printf(", nF = %3d", nF_t);
  printf("\n in the first correction magnet .: %9.2f kV/cm", ef_c1 * 1E-3);
  printf(", nF = %3d", nF_c1);
  printf("\n in the second correction magnet : %9.2f kV/cm", ef_c2 * 1E-3);
  printf(", nF = %3d", nF_c2);
  printf("\n in the charge analyzing dipole .: %9.2f kV/cm", efield * 1E-3);
  printf(", nF = %3d", nF_d);

  int nmax, ncasc;
  string fn, fn_old, fnroot;
  char answer;
  cout << "\n\n Give maximum main quantum number ........: ";
  cin >> nmax;
  cout << "\n Give number of cascade steps ............: ";
  cin >> ncasc;
  if (ncasc>nmax) ncasc=nmax;
  cout << "\n Hard or soft cut-off ? (h/s) ............: ";
  cin >> answer;
  if ((answer == 's') || (answer == 'S')) {
    // nF <= 0 signals soft cut-off
    nF_t  = -abs(nF_t);
    nF_c1 = -abs(nF_c1);
    nF_c2 = -abs(nF_c1);
    nF_d  = -abs(nF_d);
  }
  cout << "\n Give filename for output (*.fnl) ........: ";
  cin >> fnroot;

  fn = fnroot + ".hcin";
  ofstream fout(fn);
  ofstream ftxt(fnroot+".txt");
  write_inputdata_to_file(ftxt, fnroot, "TSR cooler", qion, mion, ecool, bpar * 1E7,
                          nmax, ncasc, nF_t);

  fout << "4\n";

  ftxt << "### cooler toriod" << endl;
  fn_old.clear();
  fn = fnroot + "_t";
  field_ionizing_magnet(qion, mion, ecool, z1_t, z2_t, dz_t, 100, 45, ef_t,
                        nF_t, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### 1st correction magnet" << endl;
  fn_old = fn;
  fn = fnroot + "_c1";
  field_ionizing_magnet(qion, mion, ecool, z_c1, z_c1, dz_c1, 100, 45, ef_c1,
                        nF_c1, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### 2nd correction magnet" << endl;
  fn_old = fn;
  fn = fnroot + "_c2";
  field_ionizing_magnet(qion, mion, ecool, z_c2, z_c2, dz_c2, 100, 45, ef_c2,
                        nF_c2, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### analyzing dipole magnet" << endl;
  fn_old = fn;
  fn = fnroot;
  field_ionizing_magnet(qion, mion, ecool, z_d, z_d, 0.0, rho, phi, -efield, nF_d,
                        nmax, ncasc, fn, fn_old, fout, ftxt);

  fout << "0\n";
  fout << "0\n";
  fout.close();
  ftxt.close();

  cout << "\n\n run batch job by issuing the command: nice hydrocal <" << fnroot
       << ".hcin >" << fnroot << ".log &\n";
  cout << "--------------------------------------------------------------------"
          "------------------------------\n\n";
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at TSR target
 *
 * Interactive generation of a hydrocal batch-file for the calculation of
 * field-ionization survival probabilities behind the TSR electron target
 */
void setup_batch_tsrt(void) {
  const double melectron = hydroconst::mec2_eV;

  const double rho = 115.0;  // bending radius of dipole magnet;
  const double phi = 45.0;   // bending angle of dipole magnet;
  const double z1_t = 156.8; // distance from beginning of target to toroid (in
                             // this case 50 % of field)
  const double z2_t =
      19.2; // end of target to beginning toroid (in this case 50 % of field)
  const double dz_t = 43.0; // distance within toroid (50% to 50%)
  const double z_c1 = 50.4; // distance from beginning toroid to beginning KDY2
  //  (Yes, KDY2 is nearer to the toroide than KDY1)
  const double dz_c1 = 21.0; // distance within KDY2 (2cm each side linear
                             // fringe assumed so +1cm for each side)
  const double z_c2 = 35.5;  // distance from beginning KDY2 to beginning KDY1
  const double dz_c2 = 12.0; // distance within KDY1
  const double z_d = 298.1;  // distance from beginning KDY1 to beginning dipole

  double qion, mion, ecool, bpar;

  printf("\n");
  printf("\n Setup of a batch file for the calculation of nl-specific "
         "detection probabilities at the TSR electron target");
  printf("\n");
  printf("\n Give ion charge .........................: ");
  scanf("%lf", &qion);
  printf("\n Give ion mass in u ......................: ");
  scanf("%lf", &mion);
  printf("\n Give cooling energy in eV ...............: ");
  scanf("%lf", &ecool);
  printf("\n Give magnetic guiding field in mT .......: ");
  scanf("%lf", &bpar);

  bpar *= 1E-7; // conversion of mT to Vs/(cm)^2
  double gamma = 1.0 + ecool / melectron;
  double vion = clight * sqrt(1.0 - 1.0 / (gamma * gamma)); // in cm/s
  double efield = bpar * vion;                              // in V/cm

  double ef_t = 0.5 * efield;  // EWS
  double ef_c1 = 1.9 * efield; // EWS
  double ef_c2 = 1.0 * efield; // EWS
  efield = 1.036427E-12 * gamma * vion * vion * mion / qion / rho;

  // calculate the length of the recombined ion's path through the dipole magnet
  double dz_d = rho * qion / (qion - 1.0) *
                asin((qion - 1.0) / qion * cos(phi * pi / 180.0));

  printf("\n The relativistic gamma factor is: %9.8f\n", gamma);
  printf("\n The ion velocity is ............: %9.5g cm/s\n", vion);

  printf("\n The flight times are (in the ion's frame of reference)");

  printf("\n from the beginning of the target to the toroid : %6.2f ns",
         z1_t * 1e9 / vion / gamma);
  printf("\n from the end of the target to the toroid ......: %6.2f ns",
         z2_t * 1e9 / vion / gamma);
  printf("\n from the target to the first correction magnet.: %6.2f ns",
         z_c1 * 1e9 / vion / gamma);
  printf("\n from the first to the second correction magnet.: %6.2f ns",
         z_c2 * 1e9 / vion / gamma);
  printf("\n from the second correction magnet to the dipole: %6.2f ns",
         z_d * 1e9 / vion / gamma);
  printf("\n through the target ............................: %6.2f ns",
         (z1_t - z2_t) * 1e9 / vion / gamma);
  printf("\n through the toroid ............................: %6.2f ns",
         dz_t * 1e9 / vion / gamma);
  printf("\n through the first correction magnet............: %6.2f ns",
         dz_c1 * 1e9 / vion / gamma);
  printf("\n through the second correction magnet...........: %6.2f ns",
         dz_c2 * 1e9 / vion / gamma);
  printf("\n through the dipole ............................: %6.2f ns",
         dz_d * 1e9 / vion / gamma);

  double zsum = 0.5 * (z1_t + z2_t);
  double dt_t = zsum / vion / gamma;
  printf("\n from the target center to the toroid ..........: %6.2f ns",
         dt_t * 1e9);
  zsum += z_c1;
  double dt_c1 = zsum / vion / gamma;
  printf("\n from the target center to the 1st correction m.: %6.2f ns",
         dt_c1 * 1e9);
  zsum += z_c2;
  double dt_c2 = zsum / vion / gamma;
  printf("\n from the target center to the 2nd correction m.: %6.2f ns",
         dt_c2 * 1e9);
  zsum += z_d;
  double dt_d = zsum / vion / gamma;
  printf("\n from the target center to the dipole ..........: %6.2f ns",
         dt_d * 1e9);

  int nF_t = nf(ef_t, qion);
  int nF_c1 = nf(ef_c1, qion);
  int nF_c2 = nf(ef_c2, qion);
  int nF_d = nf(efield, qion);
  printf("\n\n The motional electric fields and cut-off quantum numbers are ");
  printf("\n in the toroid ..................: %9.2f kV/cm", ef_t * 1E-3);
  printf(", nF = %3d", nF_t);
  printf("\n in the first correction magnet .: %9.2f kV/cm", ef_c1 * 1E-3);
  printf(", nF = %3d", nF_c1);
  printf("\n in the second correction magnet : %9.2f kV/cm", ef_c2 * 1E-3);
  printf(", nF = %3d", nF_c2);
  printf("\n in the charge analyzing dipole .: %9.2f kV/cm", efield * 1E-3);
  printf(", nF = %3d", nF_d);

  int nmax, ncasc;
  string fn, fn_old, fnroot;
  char answer;
  cout << "\n\n Give maximum main quantum number ........: ";
  cin >> nmax;
  cout << "\n Give number of cascade steps ............: ";
  cin >> ncasc;
  if (ncasc>nmax) ncasc=nmax;
  cout << "\n Hard or soft cut-off ? (h/s) ............: ";
  cin >> answer;
  if ((answer == 's') || (answer == 'S')) {
    // nF <= 0 signals soft cut-off
    nF_t  = -abs(nF_t);
    nF_c1 = -abs(nF_c1);
    nF_c2 = -abs(nF_c1);
    nF_d  = -abs(nF_d);
  }
  cout << "\n Give filename for output (*.fnl) ........: ";
  cin >> fnroot;

  fn = fnroot + ".hcin";
  ofstream fout(fn);
  ofstream ftxt(fnroot+".txt");
  write_inputdata_to_file(ftxt, fnroot, "TSR target", qion, mion, ecool, bpar * 1E7,
                          nmax, ncasc, nF_t);

  fout << "4\n";

  ftxt << "### cooler toroid" << endl;
  fn_old.clear();
  fn = fnroot + "_t";
  field_ionizing_magnet(qion, mion, ecool, z1_t, z2_t, dz_t, 100, 45, ef_t,
                        nF_t, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### 1st correction magnet" << endl;
  fn_old = fn;
  fn = fnroot + "_c1";
  field_ionizing_magnet(qion, mion, ecool, z_c1, z_c1, dz_c1, 100, 45, ef_c1,
                        nF_c1, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### 2nd correction magnet" << endl;
  fn_old = fn;
  fn = fnroot + "_c2";
  field_ionizing_magnet(qion, mion, ecool, z_c2, z_c2, dz_c2, 100, 45, ef_c2,
                        nF_c2, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### analyzing dipole magnet" << endl;
  fn_old = fn;
  fn = fnroot;
  field_ionizing_magnet(qion, mion, ecool, z_d, z_d, 0.0, rho, phi, -efield, nF_d,
                        nmax, ncasc, fn, fn_old, fout, ftxt);

  fout << "0\n";
  fout << "0\n";
  fout.close();
  ftxt.close();

  cout << "\n\n run batch job by issuing the command: nice hydrocal <" << fnroot
       << ".hcin >" << fnroot << ".log &\n";
  cout << "--------------------------------------------------------------------"
          "------------------------------\n\n";
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at ESR cooler
 *
 * Interactive generation of a hydrocal batch-file for the calculation of
 * field-ionization survival probabilities behind the ESR electron cooler
 */
void setup_batch_esr(void) {
  const double melectron = hydroconst::mec2_eV;

  const double rho = 625.0;    // bending radius of dipole magnet;
  const double phi = 60.0;     // bending angle of dipole magnet;
  const double z1_t = 295.0;   // distance from beginning of cooler to toroid
  const double z2_t = 45.0;    // distance from end of cooler to toroid
  const double dz_t = 60.0;    // distance within toroid
  const double z_d1 = 768.0;   // distance from toroid to 1st dipole
  const double z_d2 = 1443.75; // distance from 1st dipole to 2nd dipole
  const double z_d3 = 1443.75; // distance from 2st dipole to 3rd dipole
  const double z_d4 = 2530.50; // distance from 3rd dipole to 4th dipole
  const double dz_d = 654.50;  // distance within each dipole

  int ndip;
  double qion, mion, ecool, bpar;

  printf("\n");
  printf("\n Setup of a batch file for the calculation of nl-specific "
         "detection probabilities at the ESR electron cooler");
  printf("\n");
  printf("\n Give ion charge ..............................: ");
  scanf("%lf", &qion);
  printf("\n Give ion mass in u ...........................: ");
  scanf("%lf", &mion);
  printf("\n Give cooling energy in eV ....................: ");
  scanf("%lf", &ecool);
  printf("\n Give magnetic guiding field in mT ............: ");
  scanf("%lf", &bpar);
  printf("\n Give number of bending dipoles before detector: ");
  scanf("%d", &ndip);

  bpar *= 1E-7; // conversion of mT to Vs/(cm)^2
  double gamma = 1.0 + ecool / melectron;
  double vion = clight * sqrt(1.0 - 1.0 / pow(gamma, 2)); // in cm/s
  double efield = bpar * vion;                            // in V/cm

  double ef_t = 0.5 * efield; // 0.36*efield
  efield = 1.036427E-12 * gamma * vion * vion * mion / qion / rho;

  // calculate the length of the recombined ion's path through the dipole magnet
  // double  dz_d = rho*z/(z-1.0)*asin((z-1.0)/z*cos(phi*pi/180.0));

  printf("\n The relativistic gamma factor is: %9.8f\n", gamma);
  printf("\n The ion velocity is ............: %9.5g cm/s\n", vion);

  printf("\n The flight times are (in the ion's frame of reference)");

  printf("\n from the beginning of the cooler to the toroid : %6.2f ns",
         z1_t * 1e9 / vion / gamma);
  printf("\n from the end of the cooler to the toroid ......: %6.2f ns",
         z2_t * 1e9 / vion / gamma);
  printf("\n from the toroid to the 1st dipole .............: %6.2f ns",
         z_d1 * 1e9 / vion / gamma);
  if (ndip > 1) {
    printf("\n from the 1st to the 2nd dipole ................: %6.2f ns",
           z_d2 * 1e9 / vion / gamma);
  }
  if (ndip > 2) {
    printf("\n from the 2nd to the 3rd dipole ................: %6.2f ns",
           z_d3 * 1e9 / vion / gamma);
  }
  if (ndip > 3) {
    printf("\n from the 3rd to the 4th dipole ................: %6.2f ns",
           z_d4 * 1e9 / vion / gamma);
  }
  printf("\n through the cooler ............................: %6.2f ns",
         (z1_t - z2_t) * 1e9 / vion / gamma);
  printf("\n through the toroid ............................: %6.2f ns",
         dz_t * 1e9 / vion / gamma);
  printf("\n through the dipole ............................: %6.2f ns",
         dz_d * 1e9 / vion / gamma);

  double zsum = 0.5 * (z1_t + z2_t);
  double dt_t = zsum / vion / gamma;
  printf("\n from the cooler center to the toroid ..........: %6.2f ns",
         dt_t * 1e9);
  zsum += z_d1;
  double dt_d1 = zsum / vion / gamma;
  printf("\n from the cooler center to the 1st dipole ......: %6.2f ns",
         dt_d1 * 1e9);
  double dt_d2, dt_d3, dt_d4;
  if (ndip > 1) {
    zsum += z_d2;
    dt_d2 = zsum / vion / gamma;
    printf("\n from the cooler center to the 2nd dipole ......: %6.2f ns",
           dt_d2 * 1e9);
  }
  if (ndip > 2) {
    zsum += z_d3;
    dt_d3 = zsum / vion / gamma;
    printf("\n from the cooler center to the 3rd dipole ......: %6.2f ns",
           dt_d3 * 1e9);
  }
  if (ndip > 3) {
    zsum += z_d4;
    dt_d4 = zsum / vion / gamma;
    printf("\n from the cooler center to the 4th dipole ......: %6.2f ns",
           dt_d4 * 1e9);
  }

  int nF_t = nf(ef_t, qion);
  int nF_d = nf(efield, qion);
  printf("\n\n The motional electric fields and cut-off quantum numbers are ");
  printf("\n in the toroid ..................: %9.2f kV/cm", ef_t * 1E-3);
  printf(", nF = %3d", nF_t);
  printf("\n in the charge analyzing dipole .: %9.2f kV/cm", efield * 1E-3);
  printf(", nF = %3d", nF_d);

  char answer;
  cout << "\n\n Consider toroid ? (y/n) .................: ";
  cin >> answer;
  bool toroid_flag = ((answer == 'y') || (answer == 'Y'));

  int nmax, ncasc;
  string fn, fn_old, fnroot;
  ofstream fout;
  cout << "\n\n Give maximum main quantum number ........: ";
  cin >> nmax;
  cout << "\n Give number of cascade steps ............: ";
  cin >> ncasc;
  if (ncasc>nmax) ncasc=nmax;
  cout << "\n Hard or soft cut-off ? (h/s) ............: ";
  cin >> answer;
  if ((answer == 's') || (answer == 'S')) {
    // nF <= 0 signals soft cut-off
    nF_t  = -abs(nF_t);
    nF_d  = -abs(nF_d);
  }
  cout << "\n Give filename for output (*.fnl) ........: ";
  cin >> fnroot;
  fn = fnroot + ".hcin";
  fout = ofstream(fn);
  ofstream ftxt(fnroot+".txt");
  write_inputdata_to_file(ftxt, fnroot, "ESR", qion, mion, ecool, bpar * 1E7, nmax,
                          ncasc, nF_t, ndip);

  fout << "4\n";

  double z1, z2;
  fn_old.clear();
  if (toroid_flag) {
    ftxt << "### cooler toroid demerger" << endl;
    fn = fnroot + "_t";
    field_ionizing_magnet(qion, mion, ecool, z1_t, z2_t, dz_t, 100, 45, ef_t,
                          nF_t, nmax, ncasc, fn, fn_old, fout, ftxt);
    fn_old = fn;
    z1 = 0.0;
    z2 = 0.0;
  } else {
    z1 = z1_t;
    z2 = z2_t;
  }

  fn = fnroot;
  if (ndip > 1) {
    ftxt << "### 1st dipole magnet" << endl;
    fn += "_d1";
  } else {
    ftxt << "### analyzing dipole magnet" << endl;
    fn += "_d";
  }
  field_ionizing_magnet(qion, mion, ecool, z_d1 + z1, z_d1 + z2, dz_d, rho, phi,
                        -efield, nF_d, nmax, ncasc, fn, fn_old, fout, ftxt);

  if (ndip < 2) 
    goto finish;

  if (ndip < 3) {
    ftxt << "### analyzing dipole magnet" << endl;
  }
  else {
    ftxt << "### 2nd dipole magnet" << endl;
  }
  fn_old = fn;
  fn = fnroot + "_d2";
  field_ionizing_magnet(qion, mion, ecool, z_d2, z_d2, dz_d, rho, phi, -efield,
                        nF_d, nmax, ncasc, fn, fn_old, fout, ftxt);
  if (ndip < 3) 
    goto finish;


  if (ndip < 4) {
    ftxt << "### analyzing dipole magnet" << endl;
  }
  else {
    ftxt << "### 3rd dipole magnet" << endl;
  }
  fn_old = fn;
  fn = fnroot + "_d3";
  field_ionizing_magnet(qion, mion, ecool, z_d3, z_d3, dz_d, rho, phi, -efield,
                        nF_d, nmax, ncasc, fn, fn_old, fout, ftxt);

  if (ndip < 4)
    goto finish;

  ftxt << "### analyzing dipole magnet" << endl;
  fn_old = fn;
  fn = fnroot + "_d4";
  field_ionizing_magnet(qion, mion, ecool, z_d4, z_d4, dz_d, rho, phi, -efield,
                        nF_d, nmax, ncasc, fn, fn_old, fout, ftxt);

finish:

  fout << "0\n";
  fout << "0\n";
  fout.close();
  ftxt.close();

  cout << "\n\n run batch job by issuing the command: nice hydrocal <" << fnroot
       << ".hcin >" << fnroot << ".log &\n";
  cout << "--------------------------------------------------------------------"
          "------------------------------\n\n";
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at CRYRING cooler
 *
 * Interactive generation of a hydrocal batch-file for the calculation of
 * field-ionization survival probabilities behind the CRYRING electron cooler
 *
 * @param GSI_flag if true, length changes introduced at GSI are accounted for
 */
void setup_batch_cryring(bool GSI_flag = false) {
  const double melectron = hydroconst::mec2_eV;

  const double rho = 120.0;  // bending radius of dipole magnet;
  const double phi = 30.0;   // bending angle of dipole magnet;
  const double z1_t = 131.0; // distance from beginning of cooler to toroid
  const double z2_t = 27.0;  // distance from end of cooler to toroid
  const double dz_t = 39.0;  // distance within toroid
  const double z_c = 59.0;   // distance from toroid to correction magent
  const double dz_c = 21.0;  // distance within correction magnet
  double z_d = 47.0;         // distance from correction magnet to dipole
  if (GSI_flag)
    z_d += 21;

  double qion, mion, ecool, bpar;

  printf("\n");
  printf("\n Setup of a batch file for the calculation of nl-specific "
         "detection probabilities at the cryring electron "
         "cooler");
  printf("\n");
  printf("\n Give ion charge .........................: ");
  scanf("%lf", &qion);
  printf("\n Give ion mass in u ......................: ");
  scanf("%lf", &mion);
  printf("\n Give cooling energy in eV ...............: ");
  scanf("%lf", &ecool);
  printf("\n Give magnetic guiding field in mT .......: ");
  scanf("%lf", &bpar);

  bpar *= 1E-7; // conversion of mT to Vs/(cm)^2
  double gamma = 1.0 + ecool / melectron;
  double vion = clight * sqrt(1.0 - 1.0 / (gamma * gamma)); // in cm/s
  double efield = bpar * vion;                              // in V/cm

  double ef_t = 0.5 * efield;
  double ef_c =
      GSI_flag
          ? 1.28 * efield
          : 1.415 *
                efield; // factor 1.28 confirmed by Claude Krantz on 2026-03-12
  efield = 1.036427E-12 * gamma * vion * vion * mion / qion / rho;

  // calculate the length of the recombined ion's path through the dipole magnet
  double dz_d = rho * qion / (qion - 1.0) *
                asin((qion - 1.0) / qion * cos(phi * pi / 180.0));

  printf("\n The relativistic gamma factor is: %9.8f\n", gamma);
  printf("\n The ion velocity is ............: %9.5g cm/s\n", vion);

  printf("\n The flight times are (in the ion's frame of reference)");

  printf("\n from the beginning of the cooler to the toroid : %6.2f ns",
         z1_t * 1e9 / vion / gamma);
  printf("\n from the end of the cooler to the toroid ......: %6.2f ns",
         z2_t * 1e9 / vion / gamma);
  printf("\n from the toroid to the correction magnet ......: %6.2f ns",
         z_c * 1e9 / vion / gamma);
  printf("\n from the correction magnet to the dipole ......: %6.2f ns",
         z_d * 1e9 / vion / gamma);
  printf("\n through the cooler ............................: %6.2f ns",
         (z1_t - z2_t) * 1e9 / vion / gamma);
  printf("\n through the toroid ............................: %6.2f ns",
         dz_t * 1e9 / vion / gamma);
  printf("\n through the correction magnet .................: %6.2f ns",
         dz_c * 1e9 / vion / gamma);
  printf("\n through the dipole ............................: %6.2f ns",
         dz_d * 1e9 / vion / gamma);

  double zsum = 0.5 * (z1_t + z2_t);
  double dt_t = zsum / vion / gamma;
  printf("\n from the cooler center to the toroid ..........: %6.2f ns",
         dt_t * 1e9);
  zsum += z_c;
  double dt_c = zsum / vion / gamma;
  printf("\n from the cooler center to the correction magn. : %6.2f ns",
         dt_c * 1e9);
  zsum += z_d;
  double dt_d = zsum / vion / gamma;
  printf("\n from the cooler center to the dipole ..........: %6.2f ns",
         dt_d * 1e9);

  int nF_t = nf(ef_t, qion);
  int nF_c = nf(ef_c, qion);
  int nF_d = nf(efield, qion);
  printf("\n\n The motional electric fields are ");
  printf("\n in the toroid ..................: %9.2f kV/cm", ef_t * 1E-3);
  printf(", nF = %3d", nF_t);
  printf("\n in the correction magnet .......: %9.2f kV/cm", ef_c * 1E-3);
  printf(", nF = %3d", nF_c);
  printf("\n in the charge analyzing dipole .: %9.2f kV/cm", efield * 1E-3);
  printf(", nF = %3d", nF_d);

  int nmax, ncasc;
  string fn, fn_old, fnroot;
  char answer;

  cout << "\n\n Give maximum main quantum number ........: ";
  cin >> nmax;
  cout << "\n Give number of cascade steps ............: ";
  cin >> ncasc;
  if (ncasc>nmax) ncasc=nmax;
  cout << "\n Hard or soft cut-off ? (h/s) ............: ";
  cin >> answer;
  if ((answer == 's') || (answer == 'S')) {
    // nF <= 0 signals soft cut-off
    nF_t = -abs(nF_t);
    nF_c = -abs(nF_c);
    nF_d = -abs(nF_d);
  }
  cout << "\n Give filename for output (*.fnl) ........: ";
  cin >> fnroot;
  fn = fnroot + ".hcin";
  ofstream fout(fn);
  string storage_ring = "CRYRING";
  if (GSI_flag)
    storage_ring += "@ESR";
  ofstream ftxt(fnroot+".txt");
  write_inputdata_to_file(ftxt, fnroot, storage_ring, qion, mion, ecool, bpar * 1E7,
                          nmax, ncasc, nF_t);

  fout << "4\n";

  ftxt << "### cooler toroid demerger" << endl;
  fn_old.clear();
  fn = fnroot + "_t";
  field_ionizing_magnet(qion, mion, ecool, z1_t, z2_t, dz_t, 100, 45, ef_t,
                        nF_t, nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### correction magnet" << endl;
  fn_old = fn;
  fn = fnroot + "_c";
  field_ionizing_magnet(qion, mion, ecool, z_c, z_c, dz_c, 100, 45, ef_c, nF_c,
                        nmax, ncasc, fn, fn_old, fout, ftxt);

  ftxt << "### analyzing dipole magnet" << endl;
  fn_old = fn;
  fn = fnroot;
  field_ionizing_magnet(qion, mion, ecool, z_d, z_d, 0.0, rho, phi, -efield, nF_d,
                        nmax, ncasc, fn, fn_old, fout, ftxt);

  fout << "0\n";
  fout << "0\n";
  fout.close();
  ftxt.close();

  cout << "\n\n run batch job by issuing the command: nice hydrocal <" << fnroot
       << ".hcin >" << fnroot << ".log &\n";
  cout << "--------------------------------------------------------------------"
          "------------------------------\n\n";
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Entry to interactive batch-file setup for field-ionization
 * calculations
 */
void setup_batch(void) {
  int choice;

start:
  cout << "\n Which storage ring (target device)?\n";
  cout << "\n 1) TSR cooler";
  cout << "\n 2) TSR target";
  cout << "\n 3) ESR";
  cout << "\n 4) CRYRING";
  cout << "\n 5) CRYRING@ESR";
  cout << "\n 6) CSR";
  cout << "\n 7) CSRm\n";
  cout << "\n Make your choice : ";
  cin >> choice;

  switch (choice) {
  case 1:
    setup_batch_tsrc();
    break;
  case 2:
    setup_batch_tsrt();
    break;
  case 3:
    setup_batch_esr();
    break;
  case 4:
    setup_batch_cryring(false);
    break;
  case 5:
    setup_batch_cryring(true);
    break;
  case 6:
    setup_batch_csr();
    break;
  case 7:
    setup_batch_csrm();
    break;
  default:
    goto start;
  }
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Reads previously calculated survival probabilities from file
 */
int readfraction(vector<double> &fraction, string &header, int nmax) {
  string fn;
  ifstream fin;

  do {
    cout << "\n Give filename for fraction table (*.fnl) ..: ";
    cin >> fn;
    fn += ".fnl";
    fin.open(fn);
  } while (!fin);
  header = fn;
  header += ": ";
  string dummy;
  getline(fin, dummy);
  header += dummy;
  getline(fin, dummy);
  getline(fin, dummy);
  int n = 0, l;
  double val1, val2, val3;
  while (n <= nmax) {
    if (!(fin >> n >> l >> val1 >> val2 >> val3))
      break;
    if (n > nmax)
      break;
    fraction[(n - 1) * n / 2 + l] *= val1;
  }
  fin.close();
  return n >= nmax ? nmax : n - 1;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Interactive calculation of field-ionization cut-off quantum numbers
 */
void ngamma(void) {
  const double melectron = hydroconst::mec2_eV;
  double z, mass, ecool, cool_dip, rho, gamma, vel, efield;
  int nF, selection;
  char answer;

  printf("\n Calculation of maximum main quantum numbers which can");
  printf("\n contribute to RR in a storage ring measurement.\n");
  printf("\n The cut off due to field ionization in the dipole magnet");
  printf("\n is calulated as  nf = [z^3/(9F)]^1/4 where F is the");
  printf("\n motional electric field in V/cm.\n");
  printf("\n Also calculated is ngamma, i.e. the maximum main quantum");
  printf("\n number of Rydberg levels which on the way from the cooler");
  printf("\n to the dipole magnet radiatively decay to below nf and");
  printf("\n therefore can contribute to RR. It is calculated by");
  printf("\n iteratively solving the equation");
  printf("\n ngamma = {0.357*z^4*g0*t*[0.481 + ln(ngamma-1)]");
  printf("*(ngamma-1)^-1.38}^1/3");
  printf("\n where g0 = 2.142e10/s, and z and t");
  printf("\n denote the effective ionic charge and the flight time");
  printf("\n between cooler and magnet, respectively");
  printf("\n (see Habiltation thesis of A. Wolf).\n\n\n");
  printf("\n Give effective ionic charge z .................: ");
  scanf("%lf", &z);
  printf("\n Give atomic mass (amu) ........................: ");
  scanf("%lf", &mass);
  printf("\n Give cooling energy (eV) ......................: ");
  scanf("%lf", &ecool);
  printf("\n Give the distance D between the cooler center and");
  printf("\n entrance of the dipole magnet and the bending");
  printf("\n radius R of the magnet\n");
  printf("                             manual input : 0\n");
  printf("         TSR values (D= 472 cm, R=115 cm) : 1\n");
  printf("         ESR values (D= 938 cm, R=625 cm) : 2\n");
  printf("     CRYRING values (D= 185 cm, R=120 cm) : 3\n");
  printf(" CRYRING@ESR values (D= 206 cm, R=120 cm) : 4\n");
  printf("        CSRm values (D=1100 cm, R=760 cm) : 5\n");
  printf("        CSRe values (D=1086 cm, R=600 cm) : 6\n\n");
  printf(" Make a selection ........................: ");
  scanf("%d", &selection);
  switch (selection) {
  case 1:
    cool_dip = 472;
    rho = 115;
    break;
  case 2:
    cool_dip = 938;
    rho = 625;
    break;
  case 3:
    cool_dip = 185;
    rho = 120;
    break;
  case 4:
    cool_dip = 206;
    rho = 120;
    break;
  case 5:
    cool_dip = 1100;
    rho = 760;
    break;
  case 6:
    cool_dip = 1086;
    rho = 600;
    break;
  default:
    printf("\n Give distance D and radius R in cm: ");
    scanf("%lf %lf", &cool_dip, &rho);
  }
  gamma = 1.0 + ecool / melectron;
  vel = clight * sqrt(1.0 - 1.0 / (gamma * gamma));
  printf("\n The ion velocity is                          v = %10.4g cm/s",
         vel);
  double tau = cool_dip / vel / gamma;
  printf("\n The flight time between cooler and magnet is t = %10.4g s\n", tau);
  efield = 1.036427E-12 * vel * vel * mass / z / rho;
  for (;;) {
    nF = nf(efield, z);
    printf("\n The electric field in the dipole magnet is   F = %10.4g V/cm",
           efield);
    printf("\n leading to a cutoff main quantum number     nF = %10d", nF);
    printf("\n\n Is the electric field ok? (y/n) .....: ");
    scanf(" %c", &answer);
    if (answer == 'y' || answer == 'Y') {
      break;
    }
    printf("\n Give electric field in V/cm .........: ");
    scanf("%lf", &efield);
  }

  int ng = ngamma(tau, z, nF);
  if (ng > nF) {
    printf("\n higher Ryberg levels up to              ngamma = %10d", ng);
    printf("\n can decay radiatively to below              nF = %10d", nF);
  } else {
    printf("\n ATTENTION: ngamma = %5d <= nF = %5d", ng, nF);
    printf("\n Use nmax = nF-1 = %5d as maximum quantum number", nF - 1);
    printf(" contributing to RR.");
  }
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief recursive calculation of cascade contributions to field-ionization
 * survival probabilities
 */
double cascade(int counter, int ncascstep, int n1, int l1, int n_zero,
               vector<int> &n_list, vector<int> &l_list, vector<double> &pd_list,
               vector<double> &lt_list, const RADRATE &hydro,
               const vector<double> &psurv, double &pcascdecay) {
  counter++;

  n_list[counter] = n1;
  l_list[counter] = l1;
  pd_list[counter] = hydro.pdecay(n1, l1);
  lt_list[counter] = hydro.life(n1, l1);

  int i, n2max;
  double prod, fcascade;

  double pdecay = 0.0;
  for (i = 0; i <= counter; i++) {
    prod = 1.0;
    for (int k = 0; k <= counter; k++) {
      if (k != i)
        prod *= 1.0 - lt_list[k] / lt_list[i];
    }
    pdecay += pd_list[i] / prod;
  }
  pcascdecay = pdecay; // fraction of flux going into cascades

  double fnl = 0.0;
  if (counter < ncascstep) { // full calc. for remaining casc.
    n2max = n1 - 1;
  } else {
    n2max = (n1 <= n_zero) ? n1 - 1 : n_zero;
  }
  for (int n2 = n2max; n2 > 0; n2--) {
    for (int l2 = l1 - 1; l2 <= l1 + 1; l2 += 2) {
      if ((l2 < 0) || (l2 >= n2))
        continue;
      double br = hydro.branch(n1, l1, n2, l2);
      if (br < 1.0E-15)
        continue;
      double pdcasc = 0.0;
      fcascade = 0.0;
      if (counter < ncascstep) {
        fcascade = cascade(counter, ncascstep, n2, l2, n_zero, n_list, l_list,
                           pd_list, lt_list, hydro, psurv, pdcasc);
      }
      fnl += br * ((pdecay - pdcasc) * psurv[n2 * (n2 - 1) / 2 + l2] + fcascade);
    } // end for (l2...)
  } // end (for (n2...)

#ifdef printcasc
  for (i = 0; i < counter; i++) {
    printf("%3d %3d -> ", n_list[i], l_list[i]);
  }
  printf(" %3d %3d : %8.5f\n", n1, l1, fnl);
#endif

  return fnl;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Build the hydrogenic radiative rate matrix from RADRATE data
 *
 * The rate matrix M describes the coupled ODE system for cascade
 * populations: dP/dt = M*P.  Diagonal elements are -1/tau_i (total
 * decay rates); off-diagonal elements are B(i->j)/tau_i (feeding
 * from state i into state j).
 *
 * @param hydro   RADRATE object holding lifetimes and branching ratios
 * @param nmax    maximum principal quantum number
 * @return nstates x nstates rate matrix (nstates = nmax*(nmax+1)/2)
 */
static matrix<double> buildRateMatrix(const RADRATE &hydro, int nmax) {
  int nstates = nmax * (nmax + 1) / 2;
  matrix<double> M(nstates, nstates, 0.0);

  for (int n1 = 2; n1 <= nmax; n1++) {
    for (int l1 = 0; l1 < n1; l1++) {
      int i = (n1 - 1) * n1 / 2 + l1;
      double tau = hydro.life(n1, l1);
      if (tau <= 0.0)
        continue;
      M(i, i) = -1.0 / tau;

      for (int n2 = 1; n2 < n1; n2++) {
        for (int l2 = l1 - 1; l2 <= l1 + 1; l2 += 2) {
          if ((l2 < 0) || (l2 >= n2))
            continue;
          double br = hydro.branch(n1, l1, n2, l2);
          if (br > 0.0) {
            int j = (n2 - 1) * n2 / 2 + l2;
            M(j, i) += br / tau;
          }
        }
      }
    }
  }
  return M;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Compute all nl-dependent detection probabilities using the
 * matrix exponential method
 *
 * Solves the coupled rate equations exactly (all cascade steps) by
 * computing exp(M*t) where M is the radiative rate matrix and t is
 * the flight time.  For an extended creation region (x1 != x2) the
 * position average is evaluated with the trapezoidal rule.
 *
 * @param hydro    RADRATE object
 * @param nmax     maximum principal quantum number
 * @param x1       distance cooler entrance to field ionizer (cm)
 * @param x2       distance cooler exit to field ionizer (cm)
 * @param vion     ion velocity (cm/s)
 * @param psurv    survival probabilities for all (n,l) states
 * @param fnl_out  on return, detection probabilities for all (n,l)
 */
static void cascade_matexp(const RADRATE &hydro, int nmax, double x1,
                           double x2, double vion,
                           const vector<double> &psurv,
                           vector<double> &fnl_out) {
  matrix<double> M = buildRateMatrix(hydro, nmax);

  if (fabs(x1 - x2) < 1e-6) {
    matrix<double> expM = matrixExp(M * (x1 / vion));
    fnl_out = expM.transpose() * psurv;
  } else {
    matrix<double> expM1 = matrixExp(M * (x1 / vion));
    matrix<double> expM2 = matrixExp(M * (x2 / vion));
    matrix<double> avg = (expM1 + expM2) * 0.5;
    fnl_out = avg.transpose() * psurv;
  }
}

////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Calculation of nl-dependent hydrogenic field-ionization survival
 * probabilities
 *
 *
 * @param iselect = 0: elaborate treatment iusing the formulas given in the
 * appendix of Schippers et al. ApJ 555, 1027 (2001)
 *
 * @param iselect = 1: simple method of Zong et al. JPB 31, 3729 (1998)
 *
 * Several output files will be generated.
 * The output files are *.fn and *.fnl contain l-averaged and l-resolved
 * survival probabilities, respectively. The output file *.fm1 contains
 * l-resolved survival probabilities in matrix form for import into graphics
 * software. The output file *.fm2 contains l-resolved survival probabilities
 * weighted with the corresponing cross section in matrix form for import into
 * graphics software. The output file *.fm3 contains l-resolved hydrogenic decay
 * probabilities in matrix form for import into graphics software.
 *
 * See appendix of Schippers et al. ApJ 555, 1027 (2001) for further details.
 */

void fractionTable(int iselect) {
  const double melectron = hydroconst::mec2_eV;
  const double atomic_mass_unit = hydroconst::muc2_eV;

  if (iselect) {
    printf("Calculation of detection probabilities of hydrogenic n,l Rydberg ");
    printf("states in storage ring recombination experiments taking into\n");
    printf("account radiative decay on the way from the electron cooler to ");
    printf("the dipole magnet.\n");
    printf("The simple formula of Zong et al., JPB 31, 3729 (1998) is used.\n");
  } else {
    printf("Calculation of detection probabilities of hydrogenic n,l Rydberg ");
    printf("states in storage ring recombination experiments taking into\n");
    printf("account radiative decay on the way from the electron cooler to ");
    printf(
        "the dipole magnet or electrostatic deflector where the (motional)\n");
    printf("electric field determines the survival probability. ");
    printf("Exact hydrogenic transition rates and semi-empirical hydrogenic ");
    printf("field ionization\n");
    printf("rates are used. Cascading can be taken into account.\n\n");
  }
  printf("The fractions are written to a file (*.fnl) that can be used ");
  printf("for the calculation of recombination cross sections or rates.\n");
  printf("l-averaged fractions are written to a separate file (*.fn). ");
  printf("The n,l-selective output file *.fnl consists of 5 columns.\n");
  printf(" col 1: main quantum number n\n");
  printf(" col 2: angular momentum quantum number l\n");
  printf(" col 3: fraction of n,l-population contributing to recombination\n");
  printf(" col 4: QM cross section for RR into n,l shell at 1e-6 eV\n");
  printf(" col 5: col3*col4/(maximum cross section per n)\n\n");
  printf("Three further outputfiles *.fm1, *.fm2 and *.fm3 contain col 3, ");
  printf("col 5, and decay probabilities in matrix form,\n");
  printf("respectively. These can be used to generate 3D-plots ");
  printf("e.g. by importing them into a matrix within Origin.\n\n");
  printf(
      "Note that the required input can be most conveniently generated with ");
  printf("option 5) from the hydrocal field-ionization group.\n");

  char answer;
  double qion, mion, z1, z2, dz, ecool, rho, phi;
  printf("\n Give effective nuclear charge (0 quits) .........: ");
  scanf("%lf", &qion);
  if (qion < 0.1)
    return;
  printf("\n Give atomic mass (amu) ..........................: ");
  scanf("%lf", &mion);
  printf("\n Give cooling energy (eV) ........................: ");
  scanf("%lf", &ecool);
  printf("%6.0f %6.2f %6.2f\n", qion, mion, ecool);
  double gamma = 1.0 + ecool / melectron;
  double vion = clight * sqrt(1.0 - 1.0 / (gamma * gamma));
  double eion = ecool * mion * atomic_mass_unit / melectron;
  printf("\n The ion energy   is %12.5g keV.", eion * 0.001);
  printf("\n The ion velocity is %12.5g cm/s.\n\n", vion);
  printf(" Input of geometry values");
  printf("\n  D1: distance between beginning of cooler");
  printf(" and dipole magnet / deflector entrance");
  printf("\n  D2: distance between end of cooler");
  printf(" and dipole magnet / deflector entrance");
  printf("\n   R: bending radius of dipole magnet / deflector");
  printf("\n PHI: deflection angle of dipole magnet / deflector\n\n");
  printf("\n Give D1, D2, and R in cm ......................: ");
  scanf("%lf %lf %lf", &z1, &z2, &rho);
  printf("\n Give Phi in deg (PHI<0 for el.-stat. defl.) ...: ");
  scanf("%lf", &phi);

  if (iselect == 1) { // use simple formula of Zong et al.
    double zmean = 0.5 * (z1 + z2);
    z1 = zmean;
    z2 = zmean;
    printf("\n Using D = %5.1f cm for the simple formula of Zong et al.\n",
           zmean);
  }

  // calculate the length of the recombined ion's path through the dipole magnet
  if (phi > 0) { //  magnetic deflection
    if (fabs(qion-1.0)<1.0E-9) {
      dz = 80;
    } else {
      dz = rho * qion / (qion - 1.0) *
           asin((qion - 1.0) / qion * cos(phi * pi / 180.0));
    }
    printf("\n The distance within the magnet is %5.3f cm.\n", dz);
  } else { // electrostatic deleflection
    dz = fabs(rho * phi * pi / 180.0);
    printf("\n The distance within the deflector is %5.3f cm.\n", dz);
  }
  printf("\n Do you want to change this value ? (y/n) ........: ");
  scanf(" %c", &answer);
  if ((answer == 'y') || (answer == 'Y')) {
    printf("\n Give the distance in cm .........................: ");
    scanf("%lf", &dz);
  }
  double dwelltime = dz / vion / gamma;
  double efield;
  if (phi > 0) { // magnetic deflection
    efield = 1.036427E-12 * gamma * vion * vion * mion / qion / rho;
  } else {                       // electrostatic deflection
    efield = eion * 0.01 / qion; // for CSR 6-deg deflector
  }
  printf("\n The (motional) electric field strength is %5.2f kV/cm.",
         efield * 0.001);
  printf("\n Do you want to change this value ? (y/n) ........: ");
  scanf(" %c", &answer);
  if ((answer == 'y') || (answer == 'Y')) {
    printf("\n Give the (motional) E-field in V/cm .............: ");
    scanf("%lf", &efield);
  }
  printf("\n The dwell time in the field is %12.5g s\n", dwelltime);
  int nF = nf(efield, qion);
  printf("\n Approximate cut-off quantum number nF =  %3d\n", nF);
  printf("\n There are two models for the survival probabilities Ps,");
  printf("\n 1. a hard cut-off, i.e., Ps=1 for n<=ncut, Ps=0 for n>ncut");
  printf("\n 2. a soft cut-off calculated with FI rates ");
  printf("by Damburg and Kolosov\n");
  printf("\n Give nF (give 0 for soft cut-off) ..............:  ");
  scanf("%d", &nF);
  bool softcut_flag = nF > 0 ? false : true;

  int nmax;
  answer = 'n';
  while (answer != 'y' && answer != 'Y') {
    printf("\n Give nmax > nF .................................: ");
    scanf("%d", &nmax);

    if (nmax > nmaxfactorial / 2.0) { // nmaxfactorial is defined in hydromath.h
      nmax = int(nmaxfactorial / 2.0);
      printf("\n nmax exceeds nmaxfactorial/2");
      printf(" which is defined in hydromath.h");
      printf("\n The calculation of all necessary Clebsch-Gordan");
      printf("\n coefficients is not prossible.");
      printf("\n nmax set to nmaxfactorial/2 = %4d.\n", nmax);
    }
    double sigma_nmax = 0.0, sigma_sum = 0.0;
    for (int n = 1; n <= nmax; n++) {
      sigma_nmax = 0.0;
      for (int l = 0; l < n; l++) {
        sigma_nmax += sigmarrqm(1e-6, qion, n, l);
      }
      sigma_sum += sigma_nmax;
    }
    printf("\n Near zero energy sigma(nmax) is %8.2g%% of",
           100 * sigma_nmax / sigma_sum);
    printf(" the total RR cross section\n");
    printf("\n Is nmax o.k. ? (y/n) ............................: ");
    scanf(" %c", &answer);
  }

  // calculation of n,l-specific survival probabilities

  int n_one, n_zero;
  int psurvdim = nmax * (nmax + 1) / 2;
  vector<double> psurv(psurvdim, 0.0);
  if (softcut_flag) {
    calc_survival(dwelltime / gamma, efield, qion, nmax, psurv, n_one, n_zero);
  } else {
    for (int n = 1; n <= nmax; n++) {
      for (int l = 0; l < n; l++)
        psurv[n * (n - 1) / 2 + l] = n <= nF ? 1 : 0;
    }
    n_one = nF;
    n_zero = nF + 1;
  }

  // survival probabilities are 1 for n<=n_one and 0 for n>=n_zero
  // explicit calculations are only required for n_one < n < n_zero.

  printf("\n Now calculating hydrogenic transition rates and decay "
         "probabilities  \n");

  RADRATE hydro(nmax, qion, vion, z1 / gamma, z2 / gamma);

  vector<double> sigma(nmax, 0.0);
  vector<double> pd_list(nmax, 0.0);
  vector<double> lt_list(nmax, 0.0);
  vector<int> n_list(nmax, 0);
  vector<int> l_list(nmax, 0);

  string filenameroot, pfn, filename;
  ifstream fin;
  ofstream fout, fout2, fmatrix1, fmatrix2, fmatrix3;

new_cascade:
  fin.close();
  fout.close();
  fout2.close();
  fmatrix1.close();
  fmatrix2.close();
  fmatrix3.close();

  int ncascstep = 0;
  bool matexp_flag = false;
  cout << "\n\n Give number of cascade steps"
       << " (max if <-1, -1 for matrix exponential) ..: ";
  cin >> ncascstep;
  if (ncascstep == -1) {
    matexp_flag = true;
  } else if (ncascstep < -1) {
    ncascstep = nmax;
  }

  bool read_previous_flag = false;
  cout << "\n Multiply with values from a previous calculation?: ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    read_previous_flag = true;

    cout << "\n Give filename of previous calculation (*.fnl) ...: ";
    cin >> pfn;
    pfn += ".fnl";

    string header;

    fin.open(pfn);
    getline(fin, header);
    cout << "\n Header of previous calculation:\n " << header << "\n";
    // overread next two lines
    getline(fin, header);
    getline(fin, header);
  }

  cout << "\n Give filename for output (*.fnl, *.fm<i> i=1-3) .: ";
  cin >> filenameroot;

  filename = filenameroot + ".fnl";
  fout.open(filename);
  cout << "\n Will write l-resolved results to " << filename << ".";
  fout << "qion=" << setw(2) << qion << ", ecool=" << fixed << setprecision(2)
       << setw(7) << ecool << ", z1=" << setw(7) << fixed << setprecision(2)
       << z1 << " cm, z2=" << setw(7) << fixed << setprecision(2) << z2
       << " cm, phi=" << setw(7) << fixed << setprecision(2) << phi;
  if (matexp_flag) {  // use matrix exponentiation to compute the full cascade
    fout << " deg, ncascstep=matexp\n";
  } else {
    fout << " deg, ncascstep=" << setw(2) << ncascstep << "\n";
  }
  if (read_previous_flag) {
    fout << " previous calculation: " << pfn;
  }
  fout << "\n   n    l            f  sigma(0 eV)      f*sigma\n";

  filename = filenameroot + ".fn";
  fout2.open(filename);
  cout << "\n Will write l-averaged survival probabilites to " << filename
       << ".";

  filename = filenameroot + ".fm1";
  fmatrix1.open(filename);
  cout << "\n Will write l-resolved survival probabilities to " << filename
       << ".";

  filename = filenameroot + ".fm2";
  fmatrix2.open(filename);
  cout << "\n Will write l-resolved survival probabilities weighted with RR "
          "cross section at 0 eV to "
       << filename << ".";

  filename = filenameroot + ".fm3";
  fmatrix3.open(filename);
  cout << "\n Will write l-resolved decay probabilities to " << filename << ".";

  cout << "\n Now calculating detection probabilities \n";

  int nstates = nmax * (nmax + 1) / 2;
  vector<double> fnl_all(nstates, 0.0);
  if (matexp_flag) {
    cascade_matexp(hydro, nmax, z1 / gamma, z2 / gamma, vion, psurv, fnl_all);
  }

  int n1, l1;

  for (n1 = 1; n1 <= nmax; n1++) {
    double fn = 0.0;
    int n12 = (n1 - 1) * n1 / 2;

    double sigma_max = 0.0;
    for (l1 = 0; l1 < n1; l1++) {
      sigma[l1] = sigmarrqm(1e-12, qion, n1, l1);
      if (sigma[l1] > sigma_max) {
        sigma_max = sigma[l1];
      }
    }

    for (l1 = 0; l1 < n1; l1++) {
      double tau = hydro.life(n1, l1);
      double pdecay = hydro.pdecay(n1, l1);
      double fnl;

      if (matexp_flag) {
        fnl = fnl_all[n12 + l1];
      } else {
        fnl = (1.0 - pdecay) * psurv[n12 + l1];

        if (((n1 == 1) || (n1 == 2)) && (l1 == 0)) {
          fnl = 1.0; // 1s and 2s states are assumed to be always detected
        } else {
          pd_list[0] = pdecay;
          lt_list[0] = tau;
          n_list[0] = n1;
          l_list[0] = l1;

          int n2max;
          if (ncascstep) {
            n2max = n1 - 1;
          } else {
            n2max = (n1 <= n_zero) ? n1 - 1 : n_zero;
          }
          for (int n2 = n2max; n2 > 0; n2--) {
            int n22 = (n2 - 1) * n2 / 2;
            for (int l2 = l1 - 1; l2 <= l1 + 1; l2 += 2) {
              if ((l2 < 0) || (l2 >= n2))
                continue;
              double pcascdecay = 0.0, fcascade = 0.0;
              if (ncascstep) {
                fcascade = cascade(0, ncascstep, n2, l2, n_zero, n_list, l_list,
                                   pd_list, lt_list, hydro, psurv, pcascdecay);
              }
              double br = hydro.branch(n1, l1, n2, l2);
              fnl += br * ((pdecay - pcascdecay) * psurv[n22 + l2] + fcascade);
            }
          }
        }
      }
      // read values from previous calculation

      int np, lp;
      double fnlp, sigmap, fsp;

      if (read_previous_flag) {
        if (!(fin >> np >> lp >> fnlp >> sigmap >> fsp)) {
          cout << "\n !!! ATTENTION !!!";
          cout << "\n Actual calculation not compatible with previous one!";
          cout << "\n Program terminated at n = " << setw(4) << n1 - 1 << ".";
          fin.close();
          fout.close();
          fmatrix1.close();
          fmatrix2.close();
          fmatrix3.close();
          break;
        }
        fnl *= fnlp;
      }

      fout << defaultfloat << setw(4) << n1 << " " << setw(4) << l1 << " "
           << uppercase << setprecision(4) << setw(12) << fnl << " " << setw(12)
           << uppercase << setprecision(4) << sigma[l1] << " " << setw(12)
           << uppercase << setprecision(4) << fnl * sigma[l1] / sigma_max
           << "\n";
      fmatrix1 << setprecision(5) << setw(12) << fnl;
      fmatrix2 << setprecision(5) << setw(12) << fnl * sigma[l1] / sigma_max;
      fmatrix3 << setprecision(5) << setw(12) << pdecay;
      fn += (2.0 * l1 + 1.0) * fnl;
    }

    fn /= n1 * n1;
    fout2 << setw(5) << n1 << " " << uppercase << setprecision(4) << setw(12)
          << fn << "\n";

    for (int l = n1; l < nmax; l++) {
      fmatrix1 << setprecision(5) << setw(12) << 0.0;
      fmatrix2 << setprecision(5) << setw(12) << 0.0;
      fmatrix3 << setprecision(5) << setw(12) << 0.0;
    }
    fmatrix1 << "\n";
    fmatrix2 << "\n";
    fmatrix3 << "\n";

    if ((n1 % 10) == 0) {
      cout << "|";
    } else {
      cout << ".";
    }
    cout.flush();
  }
  fin.close();
  fout.close();
  fout2.close();
  fmatrix1.close();
  fmatrix2.close();
  fmatrix3.close();

  cout << "\n\n Another number of cascade steps? (y/n) : ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y'))
    goto new_cascade;
}

//////////////////////////////////////////////////////////////////////////
