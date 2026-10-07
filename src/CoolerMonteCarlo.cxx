/**
 * @file CoolerMonteCarlo.cxx
 *
 * @brief Monte-Carlo simulation of merged-beams recombination rate coefficients
 *
 * Reads DR cross sections from file such as files containing lorentzian peak
 * parameters or binned cross sections from autostructure. The resulting
 * rate coefficients from the Monte-Carlo convolution are written to an output
file.
 *
 * It has been tested that this MonteCarlo convolution yields the same result as
 * the convolution of the analytical cooler electron-energy distribution.
 *
 * For details see Wang et al., Eur. Phys. J. D 78, 122 (2024);
 * https://doi.org/10.1140/epjd/s10053-024-00914-7 and,
 * specifically for CRYRING, Brandau et al., Chin. Phys. C 49, 64001 (2025);
 * https://doi.org/10.1088/1674-1137/adbf81.
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
#include <chrono> // clocks and time
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional> // for "bind (random_generator, random_distribution)"
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random> // random number generators and distributions
#include "stdin_guard.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Entry to Monte Carlo simulation of electron-ion merged-beams rate
 * coefficients
 */
void CoolerMonteCarlo(void) {
  const double cLight = hydroconst::clight_cm_s;
  const double Melectron = hydroconst::mec2_eV;
  const double Matomic = hydroconst::muc2_eV;

  cout << endl
       << endl
       << "**** Merged-beam rate coefficients from  Monte-Carlo convolution of "
          "theoretical cross sections ****"
       << endl;
  cout << endl
       << " The Monte-Carlo convolution is available for the following storage "
          "rings:";
  int storage_ring_id =
      COOLER::SelectStorageRing(); // see enum StorageRingIDs in Cooler.h
  COOLER cooler(storage_ring_id);  // initializes cooler dimensions

  if (!cooler.Initialized()) {
    cout << endl << "ERROR:: Invalid choice of storage ring!" << endl << endl;
    return;
  }

  fstream fout;
  char answer;
  double Ecool, ion_mass, ion_charge, ion_deltap_over_p, ion_ktperp;
  cout << endl << " Give cooling energy (eV)  ..............................: ";
  cin >> Ecool;
  cout << endl << " Give ion mass (u) ......................................: ";
  cin >> ion_mass;
  cout << endl << " Give ion charge ........................................: ";
  cin >> ion_charge;
  cout << endl << " Give relative momentum deltaP/P (FWHM) of ion beam .....: ";
  cin >> ion_deltap_over_p;
  cout << endl << " Give transverse ion temperature (meV) ..................: ";
  cin >> ion_ktperp;
  ion_ktperp *= 0.001; // meV -> eV

  cooler.SetCoolingEnergy(Ecool); // important to do this asap
  cooler.Plot();
  double ion_gamma = 1.0 + Ecool / Melectron;
  double ion_beta = sqrt(1.0 - 1.0 / (ion_gamma * ion_gamma));
  double ion_delta_beta =
      ion_deltap_over_p * ion_beta *
      (1.0 -
       ion_beta * ion_beta); // from p=m*gamma*v=m*c*beta/sqrt(1-beta*beta)

  const double mec2 = hydroconst::mec2_eV;
  const double muc2 = hydroconst::muc2_eV;
  double mic2 = ion_mass * muc2;
  double Mratio = mec2 / mic2;
  double Mratio1 = 1.0 + Mratio;
  cooler.SetIonChargeMassRatio(ion_charge / mic2);

  double ele_ktpar, ele_ktperp;
  cout << endl << " Give electron beam temperatures kTpar and kTperp (meV) .: ";
  cin >> ele_ktpar >> ele_ktperp;
  ele_ktpar *= 0.001;  // meV -> eV
  ele_ktperp *= 0.001; // meV -> eV

  // ******************************* read DR cross section from file
  // **********************************

  int nDR = 1, theo_format = THEO_FORMAT_UNDEFINED, number_of_levels = 1,
      level_number = 1;
  string theo_filename;
  bool RRonly_flag = false;

  cout << "\n Theoretical DR data should be provided as a list of Lorentzian "
          "peak";
  cout << "\n parameters eres(eV) strength(cm2 eV) width(eV) separated by "
          "spaces";
  cout << "\n or as binned cross sections from AUTOSTRUCTURE (file extension "
          "*.bdr)";
  cout << "\n Use 'none' as input if you want only a RR calculation.\n";
  while ((theo_format != THEO_FORMAT_LORENTZ) &&
         (theo_format != THEO_FORMAT_STEIH) &&
         (theo_format != THEO_FORMAT_AUTOSDR) &&
         (theo_format != THEO_FORMAT_JAC)) {
    cout << endl
         << " Give name of DR theory data file .......................: ";
    cin >> theo_filename;
    if (theo_filename == "none") {
      theo_format = THEO_FORMAT_UNDEFINED;
      RRonly_flag = true;
      break;
    }
    nDR = open_DRtheory(theo_filename, theo_format, number_of_levels);
    if (nDR < 0) {
      cout << endl << " File " << theo_filename << " not found." << endl;
      theo_format = THEO_FORMAT_UNDEFINED;
    }
  }
  vector<double> eDR, sDR, wDR;
  double eDRmin = 0.0, eDRmax = 0.0, eDRshift = 0.0;
  if (!RRonly_flag) {
    cout << endl << " " << nDR << " DR data sets found.";
    cout << endl
         << " Overall energy shift to be applied to DR theory in eV ..: ";
    cin >> eDRshift;
    if (nDR == 0)
      return;
    cout << endl
         << " Data for " << number_of_levels
         << " different initial levels found.";
    cout << endl
         << " Give the number (in the range 1 - " << number_of_levels
         << ") of the level to be used : ";
    cin >> level_number;
    if ((level_number < 1) || (level_number > number_of_levels)) {
      cout << endl << "\n ERROR: Level number out of range!" << endl << endl;
      exit(0);
    }
    eDR.resize(nDR);
    sDR.resize(nDR);
    wDR.resize(nDR);
    read_DRtheory(theo_filename, theo_format, level_number, eDR, sDR, wDR, nDR,
                  eDRmin, eDRmax, eDRshift);
    cout << endl
         << endl
         << " " << nDR << " DR cross sections found in the CM energy range "
         << eDRmin << " - " << eDRmax << " eV";
  }
  double eCMspread = sqrt(ele_ktperp * log(2.0) * ele_ktperp * log(2.0) +
                          16.0 * log(2) * ele_ktpar * eDRmin);
  cout << endl
       << endl
       << " Minmum electron energy spread in CM system: " << eCMspread << " eV"
       << endl;

  double eCMmin, eCMmax, eCMbinwidth;
  cout << endl
       << " When defining the internal energy binnig we should take into "
          "consideration that";
  cout << endl
       << "     i) the energy bin size should be at least a factor of 10 "
          "smaller than the minimum energy spread,";
  cout << endl
       << "    ii) we might need a larger internal energy range than we "
          "finally want to use for output.";
  int nbin, nbin_add_lo = 0, nbin_add_hi = 0;
  if ((!RRonly_flag) &&
      (theo_format ==
       THEO_FORMAT_AUTOSDR)) { // cross sections are already binned
    eCMbinwidth = eDR[1] - eDR[0];
    eCMmin = eDR[0];
    eCMmax = eDR[nDR - 1] + eCMbinwidth;
    cout << endl << " Binned cross sections from AUTOSTRUCTURE";
    cout << endl
         << " Widths of cross-section bins in CM frame ...: " << eCMbinwidth
         << "eV";
    cout << endl << " CM energy range: " << eCMmin << " - " << eCMmax << " eV";
    cout << endl << " Do you need to extend the energy range ? (y/n)";
    cin >> answer;
    if ((answer == 'y') ||
        (answer == 'Y')) { // add bins at low and high energies
      double eCMmin_new, eCMmax_new;
      cout << endl
           << " Give min and max CM energy in eV .......................: ";
      cin >> eCMmin_new >> eCMmax_new;
      if (eCMmin_new < eCMmin) {
        nbin_add_lo = (int)((eCMmin - eCMmin_new) / eCMbinwidth);
        eCMmin = eCMmin - nbin_add_lo * eCMbinwidth;
      }
      if (eCMmax_new > eCMmax) {
        nbin_add_hi = (int)((eCMmax_new - eCMmax) / eCMbinwidth);
        eCMmax = eCMmax + nbin_add_hi * eCMbinwidth;
      }
    }
    nbin = nDR + nbin_add_lo + nbin_add_hi;
  } else { // cross sections are provided as a list of Lorentzian peak
           // parameters
    cout << endl
         << " Give min and max CM energy in eV .......................: ";
    cin >> eCMmin >> eCMmax;
    cout << endl
         << " Give bin width in CM system in eV ......................: ";
    cin >> eCMbinwidth;
    if (eCMbinwidth > 0.0) {
      nbin = (int)((eCMmax - eCMmin) /
                   eCMbinwidth); // number of energy bins with cross sections
    } else {
      cout << "ERROR:: Binwidth <= 0.0!!" << endl;
      exit(0);
    }
  }
  if (nbin < 1) {
    cout << "ERROR:: Too few energy bins!!" << endl;
    exit(0);
  }

  vector<double> eCMbin(nbin, 0.0); // energy grid to be used in the convolution
  vector<double> sDRbin(nbin, 0.0);
  vector<double> sRRbin(nbin, 0.0);

  if ((theo_format == THEO_FORMAT_AUTOSDR) && (nDR > 1)) {
    for (int n = 0; n < nbin; n++) {
      eCMbin[n] = eCMmin + (n + 0.5) * eCMbinwidth;
      int nn = n - nbin_add_lo;
      if ((nn >= 0) && (nn < nDR)) {
        sDRbin[n] =
            sDR[nn] /
            eCMbinwidth; // convert binned strength to binned cross section
      }
    }
  } else {
    if (RRonly_flag) {
      for (int n = 0; n < nbin; n++) {
        double ecm = eCMmin + n * eCMbinwidth;
        eCMbin[n] =
            (eCMmin < 0.0) ? exp(log(10.0) * ecm) : ecm + 0.5 * eCMbinwidth;
        ;
      }
    } else {
      bin_lorentzian_peaks(nDR, eDR, sDR, wDR, nbin, eCMmin, eCMbinwidth,
                           eCMbin, sDRbin); // see readxsec.cxx
    }
  }
  // note that eCMbin[n] contains the center energy of bin n (for a linear
  // energy scale)
  if ((RRonly_flag) && (eCMmin < 0)) {
    cout << endl << " " << nbin << " energy bins created for ";
    cout << eCMmin << " \u2264 log(Ecm/eV) \u2264 " << eCMmax << " with "
         << 1.0 / eCMbinwidth << " bins per decade";
  } else {
    cout << endl
         << " " << nbin << " energy bins of width " << eCMbinwidth
         << " eV created for ";
    cout << eCMmin << " eV \u2264 Ecm \u2264 " << eCMmax << " eV ";
  }
  double eCMoutmax = eCMmax;
  if (!RRonly_flag) {
    cout << endl
         << " Give maximum energy in eV for output ...................: ";
    cin >> eCMoutmax;
  }
  // *************************************** RR cross section
  // **********************************

  int RRnmin, RRnmax, RRlmin, RRnele;
  double RRzeff = ion_charge;
  cout << endl
       << endl
       << " Give nmin, lmin, nmax for RR ...........................: ";
  cin >> RRnmin >> RRlmin >> RRnmax;
  cout << endl << " Give number of electrons in min n,l-subshell ...........: ";
  cin >> RRnele;
  vector<double> RRfraction(1, 0.0); // we only need RRfraction[0];
  RRfraction[0] = 1.0 - RRnele / (4.0 * RRlmin + 2.0);

  if (RRzeff > 0) {
    for (int n = 0; n < nbin; n++) {
      double RReCM = fabs(eCMbin[n]) < 1.0E-7 ? 0.01 * eCMbinwidth
                                              : eCMbin[n]; // avoid RReCM=0.0
      sRRbin[n] =
          sigmarrscl(RReCM, RRzeff, RRfraction, false, RRnmin, RRnmax, RRlmin) /
          RReCM; // note that sigmarrscl returns sigma*eCM
      // cout << eCMbin[n] << " " << sRRbin[n] << endl;
    }
    cout << endl << " Calculated hydrogenic RR cross section for Z=" << RRzeff;
    cout << endl
         << "            and for n,l from " << RRnmin << "," << RRlmin << "^"
         << RRnele + 1 << " to " << RRnmax << "," << RRnmax - 1 << endl;
  }
  // **************************** display lab-energy ranges
  // ************************************

  double eLAB1min, eLAB1max, eLAB2min, eLAB2max;
  eLAB1min = EscFromEcm(-eCMbin[nbin - 1] - eCMbinwidth, Ecool,
                        ion_mass); // from kinema.cxx
  eLAB1max = EscFromEcm(-eCMbin[0], Ecool, ion_mass);
  eLAB2min = EscFromEcm(eCMbin[0], Ecool, ion_mass);
  eLAB2max = EscFromEcm(eCMbin[nbin - 1] + eCMbinwidth, Ecool, ion_mass);
  cout << endl
       << endl
       << " CM energy range from DR/RR cross section ...: " << eCMbin[0]
       << " - " << eCMbin[nbin - 1] + eCMbinwidth << " eV";
  if ((RRonly_flag) && (eCMmin < 0)) {
    cout << endl
         << " Log scale, number of points per decade .....: "
         << 1.0 / eCMbinwidth;
  } else {
    cout << endl
         << " CM energy step size ........................: " << eCMbinwidth
         << " eV";
  }
  cout << endl << " Resulting LAB electron energy ranges:";
  cout << endl << "    " << eLAB1min << " - " << eLAB1max << " eV";
  cout << endl << "    " << eLAB2min << " - " << eLAB2max << " eV";
  cout << endl;
  cout << endl << " Electrons faster or slower than ions? (f/s) ...........: ";
  cin >> answer;
  bool fast_flag = ((answer == 'f') || (answer == 'F'));

  cout << endl << " Rate coefficient or energy distribution ? (r/d) .......: ";
  cin >> answer;
  bool distri_flag = ((answer == 'd') || (answer == 'D'));
  double eCMdistri = 0.0;
  if (distri_flag) {
    cout << endl
         << " Give the CM energy of distribution in eV ..............: ";
    cin >> eCMdistri;
  }

  // ****************** seed for the random generator
  // ***************************
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();

  int nMC;
  cout << endl << " Give the number of Monte-Carlo iterations .............: ";
  cin >> nMC;

  // *****************************  time stamp
  // **************************************
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  // ********************* open output file and write header
  // ************************
  string outfilename, extension;
  if (distri_flag) {
    extension = ".mcd";
  } else {
    extension = ".mcr";
  }
  cout << endl
       << " Give filename for output (" << extension
       << ") .......................: ";
  cin >> outfilename;
  outfilename += extension;

  fout.open(outfilename, fstream::out);
  fout << "####################################################################"
          "########################################"
          "###########"
       << endl;
  if (distri_flag) {
    fout << "### CM-energy distribution from a Monte-Carlo simulation of an "
            "electron-ion merged-beams experiment"
         << endl;
  } else {
    fout << "### Rate coefficients from a Monte-Carlo simulation of an "
            "electron-ion merged-beams experiment"
         << endl;
  }
  fout << "###" << endl;
  fout << "###            Publication : Eur. Phys. J. D 78, 122 (2024)" << endl;
  fout << "###                    DOI : https://doi.org/10.1140/epjd/s10053-024-00914-7" << endl;
  fout << "###" << endl;
  fout << "###  hydrocal revision     : " << HYDROCAL_REVISION
       << endl; // defined in buildinfo.h
  fout << "###               filename : " << outfilename << endl;
  fout << "###      start date & time : " << asctime(timeinfo);
  cooler.Print(fout);
  fout << "###            DR filename : " << theo_filename << endl;
  if (!RRonly_flag) {
    fout << "###           level number : " << level_number << endl;
    fout << "###   DR energy shift (eV) : " << eDRshift << endl;
  }
  if (fast_flag) {
    fout << "###              electrons : faster then ions" << endl;
  } else {
    fout << "###              electrons : slower then ions" << endl;
  }
  if (distri_flag) {
    fout << "### nominal CM energy (eV) : " << eCMdistri << endl;
    double Elab = fast_flag ? EscFromEcm(eCMdistri, Ecool, ion_mass)
                            : EscFromEcm(-eCMdistri, Ecool, ion_mass);
    fout << "###     -> LAB energy (eV) : " << Elab << endl;
  } else {
    fout << "###            RR Zeff     : " << RRzeff << endl;
    fout << "###            RR nmin     : " << RRnmin << endl;
    fout << "###            RR lmin     : " << RRlmin << endl;
    fout << "###            RR nele     : " << RRnele << endl;
    fout << "###            RR nmax     : " << RRnmax << endl;
  }
  fout << "###    cooling energy (eV) : " << Ecool << endl;
  fout << "###    ion mass (u)        : " << ion_mass << endl;
  fout << "###    ion charge          : " << ion_charge << endl;
  fout << "###    ion delta_p/p       : " << ion_deltap_over_p << endl;
  fout << "###    ion gamma           : " << ion_gamma << endl;
  fout << "###    ion beta            : " << ion_beta << endl;
  fout << "###    ion delta_beta/beta : " << ion_delta_beta / ion_beta << endl;
  fout << "###    ion kTperp (meV)    : " << ion_ktperp * 1000.0 << endl;
  fout << "###  electron kTpar  (meV) : " << ele_ktpar * 1000.0 << endl;
  fout << "###  electron kTperp (meV) : " << ele_ktperp * 1000.0 << endl;
  if ((RRonly_flag) && (eCMmin < 0.0)) {
    fout << "### log min CM energy (eV) : " << eCMmin << endl;
    fout << "### log max CM energy (eV) : " << eCMmax << endl;
  } else {
    fout << "###     min CM energy (eV) : " << eCMmin << endl;
    fout << "###     max CM energy (eV) : " << eCMmax << endl;
  }
  fout << "###     eCM bin width (eV) : " << eCMbinwidth << endl;
  fout << "### max output energy (eV) : " << eCMoutmax << endl;
  fout << "###     no. of energy bins : " << nbin << endl;
  fout << "###            random seed : " << seed << endl;
  fout << "### Monte-Carlo iterations : " << nMC << endl;
  if (!distri_flag) {
    fout << "##################################################################"
            "########################################"
            "###################################################"
         << endl;
    fout << "### LAB energy (eV) [1]   CM energy (eV) [2]   alpha DR (cm3/s) "
            "[3]   alpha RR (cm3/s) [4]   alpha RR+DR "
            "(cm3/s) [5]   Rate per 1e6 ions (Hz) [6]   norm [7] "
         << endl;
    fout << "###---------------------------------------------------------------"
            "----------------------------------------"
            "---------------------------------------------------"
         << endl;
  }

  cout << endl << " Additional diagnostic output on screen? (y/n/e) .......: ";
  cin >> answer;
  bool diagnostic_flag = ((answer == 'y') || (answer == 'Y'));
  bool exit_flag = ((answer == 'e') || (answer == 'E'));

  // number of threads for parallel processing
  int threads = omp_get_num_procs();
  threads = (threads > 4) ? 4 : 1;
  omp_set_num_threads(threads);

  if (!distri_flag) {
    int noutbin = (int)((eCMoutmax - eCMmin) /
                        eCMbinwidth); // number of energy bins for output
    cout << endl
         << endl
         << " Now using " << threads << " parallel threads for performing "
         << nMC << " Monte-Carlo iterations" << endl;
    cout << " for each of up to " << noutbin << " relative energies." << endl
         << " Remain patient! " << endl
         << endl;
    if ((RRonly_flag) && (eCMmin < 0))
      eCMoutmax = exp(log(10.0) * eCMmax); // if log energy scale
  }

  //************************ initialize random generators and random
  //distributions  *******************
  mt19937_64 randgen(
      seed); // 64bit "Mersenne Twister 19937" random number generator

  const double ln2 = hydroconst::sqrtln2 * hydroconst::sqrtln2;
  double cooler_length = cooler.MaxOverlapLength();
  uniform_real_distribution<double> dist_cooler(0.0, cooler_length);
  normal_distribution<double> dist_ele_perp(0.0, sqrt(ele_ktperp / Melectron));
  normal_distribution<double> dist_ele_par(0.0, sqrt(ele_ktpar / Melectron));
  normal_distribution<double> dist_ion_perp(
      0.0, sqrt(ion_ktperp / (ion_mass * Matomic)));
  normal_distribution<double> dist_ion_par(
      0.0, ion_delta_beta /
               sqrt(8.0 * ln2)); // with conversion from FWHM to Gaussian sigma

  // setup z-axis for output of mean Ecm vs. z
  const int nzpos = 200;
  double dzpos = cooler_length / ((double)nzpos);
  vector<double> dist_z_count(nzpos,
                              0.0); // this histrogram contains the number of
                                    // occurrences of a given z-position
  vector<double> dist_z_eCM(
      nzpos, 0.0); // this histrogram sums up the cm energies at each z-position
  vector<double> nele_z(nzpos, 0.0); // electron density at each z-position
  vector<double> dist_eCM(
      nbin,
      0.0); // this is the histogramm that contains the CM energy distribution
  vector<double> nele_eCM(nbin, 0.0); // electron density as function of eCM
  bool norm_ok_flag = false;
  int first_bin = 0, last_bin = 0;
  // *********** calculation of rate coefficients alpha(erel) by convolution of
  // sigma(eCM) with f(eCM,erel) *******
  for (int n = 0; n < nbin; n++) { // loop over center-of-mass energies
    int out_of_zrange_counter = 0.0;
    double eCM = distri_flag ? eCMdistri : eCMbin[n];
    if (eCM > (eCMoutmax + 0.51 * eCMbinwidth))
      break;
    double eLAB, eLABlo, eLABhi;
    // relativistic calculation of lab energy (space charge included if
    // electron_current > 0)
    eLABlo = EscFromEcm(-eCM, Ecool, ion_mass);
    eLABhi = EscFromEcm(eCM, Ecool, ion_mass);
    eLAB = fast_flag ? eLABhi : eLABlo;
    cooler.SetLaboratoryEnergy(eLAB); // important to do this asap

    double nele_mean = 0.0;
    for (int j = 0; j < nzpos; j++) { // z-dependent electron density
      nele_z[j] = cooler.BeamDensity((j + 0.5) * dzpos);
      nele_mean += nele_z[j];
    }
    nele_mean /= nzpos;

    // **************************  Monte-Carlo simulation of the CM energy
    // distribution **************
    int dist_eCM_count = 0;
    // an energy distribution is simulated specifically for each eCM, i.e. we
    // need to zero first
    for (int j = 0; j < nbin; j++) {
      dist_eCM[j] = 0.0;
      nele_eCM[j] = 0.0;
    }
    for (int i = 0; i < nMC; i++) {
      // Monte-Carlo loop generating energy distribution
      double zpos = dist_cooler(randgen);   // z-axis in beam direction
      int izpos = (int)floor(zpos / dzpos); // array index of position
      double ele_beta_x = dist_ele_perp(randgen);
      double ele_beta_y = dist_ele_perp(randgen);
      double ele_beta_z = dist_ele_par(randgen);

      double ion_gamma_cooler = ion_gamma;
      bool ok = cooler.BeamBeta(zpos, ion_gamma_cooler, ele_beta_x, ele_beta_y,
                                ele_beta_z);
      if (!ok) {
        out_of_zrange_counter++;
        continue; // if zpos out of range
      }
      double ion_beta_cooler =
          sqrt(1.0 - 1.0 / (ion_gamma_cooler * ion_gamma_cooler));
      double ion_beta_x = dist_ion_perp(randgen);
      double ion_beta_y = dist_ion_perp(randgen);
      double ion_beta_z = ion_beta_cooler + dist_ion_par(randgen);

      // Relativistic addition of velocities
      // small x, y and large z components may lead to numerical errors!!
      double ele_beta_sqr = ele_beta_x * ele_beta_x + ele_beta_y * ele_beta_y +
                            ele_beta_z * ele_beta_z;
      double ion_beta_sqr = ion_beta_x * ion_beta_x + ion_beta_y * ion_beta_y +
                            ion_beta_z * ion_beta_z;
      double ion_ele_prod = ion_beta_x * ele_beta_x + ion_beta_y * ele_beta_y +
                            ion_beta_z * ele_beta_z;
      double gamma_rel = (1.0 - ion_ele_prod) /
                         sqrt((1.0 - ele_beta_sqr) * (1.0 - ion_beta_sqr));
      double e_cm =
          mic2 * (sqrt(1.0 + Mratio * (Mratio + 2.0 * gamma_rel)) - Mratio1);

      // sorting into energy bin
      int bin;
      if ((RRonly_flag) && (eCMmin < 0)) {
        bin = (int)floor((log10(e_cm) - eCMmin) / eCMbinwidth);
      } else {
        bin = (int)floor((e_cm - eCMmin) / eCMbinwidth);
      }
      if ((bin >= 0) && (bin < nbin)) {
        nele_eCM[bin] += nele_z[izpos];
        dist_eCM[bin] += 1.0;
        dist_eCM_count++;
      }
      if (distri_flag) { // energy as function of position
        dist_z_count[izpos] += 1.0;
        dist_z_eCM[izpos] += e_cm;
      }
    } // end for(i ...) end of Monte-Carlo loop
    //*******************************************
    if (diagnostic_flag) {
      cout << endl
           << " Ecool: " << Ecool << " eV, Ecm: " << eCM
           << " eV, Elab: " << eLAB
           << " eV, number of histogram entries: " << dist_eCM_count
           << " number of out-of-z-range events: " << out_of_zrange_counter;
    }

    double valid_nMC =
        nMC - out_of_zrange_counter; // number of valid Monte Carlo iterations
                                     // ("out of z-range" events not considered)
    double dist_norm =
        dist_eCM_count / valid_nMC; // should be one or very close to one
    // dist_norm will be lower than one, if a significant part of the energy
    // distribution is outside of the interval [eCMmin, eCMmax] This will be the
    // case for eCM being close to either eCMmin or eCMmax. We will discar these
    // cases. We assume, that disnorm will be smaller than one for the first few
    // energy points and for the last few energy points. Starting from low
    // energies, the first few energy distributions will have dist_norm < 1, the
    // the norm will be 1 (ore close to 1) for many energeies until it drops
    // again towards the highest energies. As soon as the norm is above a
    // threshold we will set the norm_ok_flag. If afterwards the dist_norm drops
    // below the threshold again the program will stop iterating and terminate
    // normally.

    if (!distri_flag) {
      if (dist_norm < 0.98) { // we allow for 2% uncertainty
        if (norm_ok_flag) {
          break;
        } else {
          continue;
        }
      } else {
        if (first_bin == 0)
          first_bin = n;
        last_bin = n;
        norm_ok_flag = true;
      }
    }
    // at this point we have an energy distribution dist_eCM[]=f(eCM,erel)*deCM
    // that can be used for convolution (integration deCM) with the DR and RR
    // cross sections
    if (distri_flag) {
      fout << "### drift tube voltage (V) : " << cooler.DriftTubeVoltage()
           << endl;
      fout << "################################################################"
              "########################################"
              "################################################################"
              "########################################"
              "#################"
           << endl;
      fout << "### CM energy (eV) [1]   energy distribution (1/eV) [2]   "
              "electron_density (cm-3) [3]   binned DR cross "
              "section (cm2) [4]   binned RR cross section (cm2) [5]   z (cm) "
              "[6]   mean Ecm (eV) [7]   mean electron "
              "density (cm-3) [8]"
           << endl;
      fout << "###-------------------------------------------------------------"
              "----------------------------------------"
              "----------------------------------------------------------------"
              "----------------------------------------"
              "-----------------"
           << endl;
    }
    if (distri_flag) {
      // -------------- output of energy distribution (first part)
      // -------------------------------------
      for (int j = 0; j < nbin; j++) {
        if (dist_eCM[j] > 0.0)
          nele_eCM[j] /= dist_eCM[j];  // normalize electron density
        dist_eCM[j] /= dist_eCM_count; // normalize distribution such that the
                                       // integral yields 1
        double ecm = eCMbin[j];
        double ecm_mean = 0.0;
        // output of energy distribution, electron density and binned DR and RR
        // cross sections
        fout << ecm << ", " << dist_eCM[j] / eCMbinwidth << ", " << nele_eCM[j]
             << ", " << sDRbin[j] << ", " << sRRbin[j] << ", ";
        // output of <Ecm> and electron density vs. z position
        if (j < nzpos) {
          if (dist_z_count[j] > 0.0) {
            ecm_mean = dist_z_eCM[j] / dist_z_count[j];
          } else {
            ecm_mean = 0.0;
          }
          fout << (j + 0.5) * dzpos - 0.5 * cooler_length << ", " << ecm_mean
               << ", " << nele_z[j] << endl;
        } else {
          fout << ", , " << endl;
        }
      } // end for (j...)
      // -------------- output of energy distribution (second part)
      // -------------------------------------
      for (int j = nbin; (nbin < nzpos) && (j < nzpos); j++) {
        double ecm_mean = 0.0;
        if (dist_z_count[j] > 0.0) {
          ecm_mean = dist_z_eCM[j] / dist_z_count[j];
        }
        fout << " , , , , , " << (j + 0.5) * dzpos - 0.5 * cooler_length << ", "
             << ecm_mean << ", " << nele_z[j] << endl;
      } // end for (j...)
      break; // no loop for(n...) over energies required since the energy
             // distribution is evaluated only for one specified energy
    } // -------------------------------------------------------------------------------------------------

    //***************************************************************************************************
    //****************  Convolution of DR and RR cross sections with the CM
    //energy distribution *********
    double alphaDR = 0.0, alphaRR = 0.0;
    for (int j = 0; j < nbin; j++) {
      double ecm = eCMbin[j];
      double Mratio1ecm = Mratio1 + ecm / mic2;
      double gamma_rel =
          1.0 + 0.5 * (Mratio1ecm * Mratio1ecm - Mratio1 * Mratio1) /
                    Mratio; // from exact Ecm
      double v_rel = cLight * sqrt(1.0 - 1.0 / (gamma_rel * gamma_rel));
      /* // assuming constant electron density
      dist_eCM[j] /= dist_eCM_count; // normalize distribution such that the
      integral yields 1
      //integration alpha(E)=v_rel(eCM)*sigma(eCM)*f(eCM,E)*deCM, note that
      dist_eCM corresponds to f(eCM,E)*deCM if (!RRonly_flag) alphaDR +=
      v_rel*sDRbin[j]*dist_eCM[j]; alphaRR += v_rel*sRRbin[j]*dist_eCM[j];
      */
      // nele_eCM[j] /= dist_eCM[j]
      // dist_eCM[j] /= dist_eCM_count; // normalize distribution such that the
      // integral yields 1
      //  we consider the z-dependent electron density is in the integral!
      //  integration
      //  nele*alpha(E)=nele(eCM)*v_rel(eCM)*sigma(eCM)*f(eCM,E)*deCM, note that
      //  nele_eCM[j]*dist_eCM[j] = nele_eCM[j]/dist_eCM_count
      // and that nele_eCM corresponds to ne(E)*f(eCM,E)*deCM
      if (!RRonly_flag)
        alphaDR += v_rel * sDRbin[j] * nele_eCM[j] / dist_eCM_count;
      alphaRR += v_rel * sRRbin[j] * nele_eCM[j] / dist_eCM_count;
    } // end for(j ...), end of convolution loop
    alphaDR /= nele_mean;
    alphaRR /= nele_mean;
    double eLAB_gamma = 1.0 + eLAB / Melectron;
    double eLAB_beta = sqrt(1.0 - 1.0 / (eLAB_gamma * eLAB_gamma));
    double LoverC = cooler.NominalOverlapLength() / cooler.RingCircumference();
    double count_rate = 1E6 * (alphaRR + alphaDR) *
                        (1.0 - eLAB_beta * ion_beta) * LoverC *
                        nele_mean; // per 1E6 stored ions
    fout << eLAB << ",  " << eCM << ", " << alphaDR << ", " << alphaRR << ", "
         << alphaDR + alphaRR << ", " << count_rate << ", " << dist_norm
         << endl;
    if ((!diagnostic_flag) && (((n + 1) % 100) == 0)) {
      cout << " " << n + 1 << flush;
    }
    if ((!diagnostic_flag) && (((n + 1) % 1500) == 0)) {
      cout << endl;
    }
  } // end for(n ...), end of loop over energies

  // ***************************** cleanup ******************************
  fout.close();
  if (distri_flag) {
    cout << endl
         << endl
         << " Distribution data written to output file: " << outfilename
         << endl;
  } else {
    if (last_bin > 0) {
      cout << endl
           << endl
           << " Data for " << last_bin - first_bin + 1
           << " energies in the range ";
      cout << eCMbin[first_bin] << " - " << eCMbin[last_bin] << " eV" << endl
           << " written to output file " << outfilename << "." << endl;
    } else {
      cout << " ATTENTION: No output generated! Probably, the energy range "
              "needs to be extended (see comment below)."
           << endl;
    }
  }
  cout << endl
       << " Note that the "
          "norm"
          " in column 6 of the output file indicates whether";
  cout << endl
       << " the corresponding energy distribution was properly normalized.";
  cout << endl
       << " Ideally the "
          "norm"
          " is equal to 1.0. If it drops below 0.98";
  cout << endl
       << " the corresponding energy bin is discarded. This happens when the "
          "energy";
  cout << endl
       << " distribution extends to beyond the limits of the chosen energy "
          "range.";
  cout << endl
       << " The energy range should be extended if some of the discarded bins "
          "are required.";
  cout << endl;
  if (exit_flag) {
    exit(0);
  } else {
    return;
  }
}
