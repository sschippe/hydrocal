/**
 *  @file RRDRxsec.cxx
 *
 * @brief hydrogenic RR cross sections and DR peak cross sections
 *
 *  @par CREATION
 *  @author Stefan Schippers
 *  @date 1997
 *
 *  @par VERSION
 *  @verbatim
// SPDX-License-Identifier: MIT
 *  @endverbatim
*/

#include "RRDRxsec.h"
#include "CoolerFractions.h"
#include "DiracRate.h"
#include "DiracRR.h"
#include "fele.h"
#include "hydroconst.h"
#include "hydromath.h"
#include "osci.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// interpolation arrays
int _sigma_n = 0;
vector<double> _sigma_x, _sigma_y, _sigma_dy;

//////////////////////////////////////////////////////////////////////////
/**
 * @brief initializes interpolation arrays
 *
 * @param energy array of energy values
 * @param xesc array of cross section values
 * @param npts array size
 */
void init_sigma_interpolation(vector<double> energy, vector<double> xsec,
                              int npts) {
  _sigma_n = npts;
  _sigma_x.resize(npts);
  _sigma_y.resize(npts);
  _sigma_dy.resize(npts);
  for (int i = 0; i < npts; i++) {
    _sigma_x[i] = energy[i];
    _sigma_y[i] = xsec[i];
    //  cout << _sigma_x[i] << " " << _sigma_y[i] << endl;
  }
  spline(_sigma_n, _sigma_x, _sigma_y, _sigma_dy);
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief returns interpolated cross section
 *
 * @param e energy
 * @param . all other parameters are not used
 *
 * @return interpolated cross section times energy
 * @return 0 if the energy is outside the interpolation range
 *
 */
double sigma_interpolated(double e, double /*z*/, vector<double> /*fraction*/,
                          int /*use_fraction*/, int /*nmin*/, int /*nmax*/, int /*lmin*/) {
  if ((e < _sigma_x[0]) || (e > _sigma_x[_sigma_n - 1])) {
    return 0;
  } else {
    return e * splint(e, _sigma_n, _sigma_x, _sigma_y, _sigma_dy);
  }
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief convolution of a delta peak like cross section
 *
 * @param fele electron energy distribution
 * @param er resonance energy in eV
 * @param eV electron-ion collision energy in eV
 * @param kTpar longitudinal temperature of the electron beam in eV
 * @param kTpesp transversal temperature of the electron beam in eV
 *
 * @return alpha_DR(E)
 */
double deltapeak(FELE fele, double er, double eV, double ktpar, double ktperp,
                 double a) {
  double vsigma =
      a * hydroconst::clight_cm_s * sqrt(2.0 * er / hydroconst::mec2_eV);
  return (*fele)(eV, er, ktpar, ktperp) * vsigma;
}


//////////////////////////////////////////////////////////////////////////
/**
 * @brief lorentzian resonance cross section
 *
 * @param e electron-ion collision energy
 * @param er resonance energy in eV
 * @param fraction vector of length >=2 passing resonance parameters
 * @param fraction[0] peak area  in cm2 eV
 * @param fraction[1] peak width in eV 
 * @param eV electron-ion collision energy in eV
 * @param kTpar longitudinal temperature of the electron beam in eV
 * @param kTpesp transversal temperature of the electron beam in eV
 *
 * @return sigma*e if width>0
 * @return sigma*er if (width<0) and (area>0)
 * @return sigma*(er+(width/2)²/er) if (width<0) and (area<0)
 */
double sigmalorentzian(double e, double er, vector<double> fraction,
                       int /*use_fraction*/, int /*nmin*/, int /*nmax*/, int /*lmin*/) {
  double a = fraction[0];   // peak area
  double g = fraction[1];   // peak width
  double gg = 0.25 * g * g; // (width/2)^2
  double esigma =
      0.5 * fabs(a) * e / hydroconst::pi * fabs(g) / ((e - er) * (e - er) + gg);
  if (g < 0) {
    if (a < 0) {
      esigma *= (er + gg / er) / e;
    } else {
      esigma *= er / e;
    }
  }
  return esigma;
}

////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief test cross section yielding alpharr = 1.0
 */
double sigmarrtest(double e, double /*z*/, vector<double> /*fraction*/,
                   int /*use_fraction*/, int /*nmin*/, int /*nmax*/, int /*lmin*/) {
  return 1.0 / hydroconst::clight_cm_s / sqrt(2.0 / hydroconst::mec2_eV / e);
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief nonrelativistic hydrogenic RR cross section for subshell n,l at zero energy
 *
 * @param z nuclear charge
 * @param n principal quantum number of initial subshell
 * @param l orbital angular momentum quantum number
 *
 * @return RR cross section times electron-ion collision energy (E) in cm2 eV for E=0
 */
double sigmarrqm(double z, int n, int l) { // 
  const double Pi = hydroconst::pi;
  const double Alpha = hydroconst::alpha;
  const double A0 = hydroconst::a0_cm;
  const double Rydberg = hydroconst::Ryd_eV;
  double nfac = z / n / n * Rydberg;
  double sigma = nfac * nfac * (2.0 * l + 1.0) * fosciBC(n, l);
  return 2.0 * Pi * Pi * A0 * A0 * sigma / Rydberg * Alpha * Alpha * Alpha;
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief nonrelativistic quantum mechanical hydrogenic RR cross section for subshell n,l
 *
 * from RR cross section via detailed balance 
 *
 * @param eV electron-ion collision energy in eV
 * @param z nuclear charge
 * @param n principal quantum number of initial subshell
 * @param l orbital angular momentum quantum number
 *
 * @return RR cross section times electron-ion collision energy in cm2 eV
 */
double sigmarrqm(double eV, double z, int n, int l) { 
  const double Pi = hydroconst::pi;
  const double Alpha = hydroconst::alpha;
  const double A0 = hydroconst::a0_cm;
  const double Rydberg = hydroconst::Ryd_eV;
  double e = eV / Rydberg / z / z;
  double nfac = (eV + z * z / n / n * Rydberg) / z;
  double sigma = nfac * nfac * (2.0 * l + 1.0) * fosciBC(e, n, l);
  return 2.0 * Pi * Pi * A0 * A0 * sigma / Rydberg * Alpha * Alpha * Alpha;
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief nonrelativistic hydrogenic photoionization cross section for subshell n,l
 *
 * from RR cross section via detailed balance 
 *
 * @param Eph_eV photon energy in eV
 * @param z nuclear charge
 * @param nmin principal quantum number of the shell for which IP is specified
 * @param n principal quantum number of initial subshell
 * @param l orbital angular momentum quantum number
 * @param IP ionization potential of nmin-shell in eV
 *
 * @return photoionization cross section times photon energy in cm2 eV
 */
double sigmaPIqm(double Eph_eV, double z, int nmin, int n, int l, double IP) {
  const double Rydberg = hydroconst::Ryd_eV;
  double gf = 1;
  double gi = 4 * l + 2;
  double Erel = Eph_eV - Rydberg * z * z * (1.0 / (n * n) - 1.0 / (nmin * nmin)) + IP;
  return 1.022e6 * gf / gi / Eph_eV * sigmarrqm(Erel, z, n, l);
}

////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief nonrelativistic hydrogenic photoionization cross section summed over subshells
 *
 * @param Eph_eV photon energy in eV
 * @param z nuclear charge
 * @param fraction fractional populations for n,l subshells
 * @param use_fraction flags whether the above fractions should be used or not 
 * @param nmin principal quantum number of lowest shell
 * @param nmax prinicpal quantum number of highest shell
 * @param lmin orbital angular momentum quantum number of lowest subshell
 * @param IP ionization potential of nmin-shell in eV 
 *
 * @return photoionization cross section times photon energy in cm2 eV
 */
double sigmaPIqm(double Eph_eV, double z, vector<double> fraction, int use_fraction,
                 int nmin, int nmax, int lmin, double IP) {
  double lterm = sigmaPIqm(Eph_eV, z, nmin, nmin, lmin, IP);
  int n12 = (nmin - 1) * nmin / 2;
  if (use_fraction) {
    lterm *= fraction[n12 + lmin];
  } else {
    lterm *= fraction[0];
  }
  double lsum = lterm;
  for (int l = lmin + 1; l < nmin; l++) {
    lterm = sigmaPIqm(Eph_eV, z, nmin, nmin, l, IP);
    if (use_fraction) {
      lterm *= fraction[n12 + l];
    }
    lsum += lterm;
  }
  double nsum = lsum;
  for (int n = nmin + 1; n <= nmax; n++) {
    lsum = 0.0;
    n12 = (n - 1) * n / 2;
    for (int l = 0; l < n; l++) {
      lterm = sigmaPIqm(Eph_eV, z, nmin, n, l, IP);
      if (use_fraction) {
        lterm *= fraction[n12 + l];
      }
      lsum += lterm;
    }
    nsum += lsum;
  }
  return nsum;
}

////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief nonrelativistic quantum mechanical hydrogenic RR cross section summed over subshells
 *
 * @param eV electron-ion collision energy in eV
 * @param z nuclear charge
 * @param fraction fractional populations for n,l subshells
 * @param use_fraction flags whether the above fractions should be used or not 
 * @param nmin principal quantum number of lowest shell
 * @param nmax prinicpal quantum number of highest shell
 * @param lmin orbital angular momentum quantum number of lowest subshell
 *
 * @return RR cross section times electron-ion collision energy in cm2 eV
 */
double sigmarrqm(double eV, double z, vector<double> fraction, int use_fraction,
          int nmin, int nmax, int lmin) { 
  const double Pi = hydroconst::pi;
  const double Alpha = hydroconst::alpha;
  const double A0 = hydroconst::a0_cm;
  const double Rydberg = hydroconst::Ryd_eV;

  double e = eV / Rydberg / z / z;
  double nfac = (eV + z * z / nmin / nmin * Rydberg) / z;
  double lterm = (2.0 * lmin + 1.0) * fosciBC(e, nmin, lmin);
  int n12 = (nmin - 1) * nmin / 2;
  if (use_fraction) {
    lterm *= fraction[n12 + lmin];
  } else {
    lterm *= fraction[0];
  }
  double lsum = lterm;
  for (int l = lmin + 1; l < nmin; l++) {
    lterm = (2.0 * l + 1.0) * fosciBC(e, nmin, l);
    if (use_fraction) {
      lterm *= fraction[n12 + l];
    }
    lsum += lterm;
  }
  double nsum = nfac * nfac * lsum;
  for (int n = nmin + 1; n <= nmax; n++) {
    nfac = (eV + z * z / n / n * Rydberg) / z;
    lsum = 0.0;
    n12 = (n - 1) * n / 2;
    for (int l = 0; l < n; l++) {
      lterm = (2.0 * l + 1.0) * fosciBC(e, n, l);
      if (use_fraction) {
        lterm *= fraction[n12 + l];
      }
      lsum += lterm;
    }
    nsum += nfac * nfac * lsum;
  }
  return 2.0 * Pi * Pi * A0 * A0 * nsum / Rydberg * Alpha * Alpha * Alpha;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief relativistic (Dirac) hydrogenic RR cross section in dipole approximation for nonrelativistic subshell n,l
 *
 * @param eV electron-ion collision energy in eV
 * @param z nuclear charge
 * @param n principal quantum number of initial subshell
 * @param l orbital angular momentum quantum number
 *
 * @return RR cross section times electron-ion collision energy in cm2 eV
 */
double sigmarrdid(double eV, double z, int n, int l)
{  
  double w1 = 2.0*(l+0.5)+1.0;
  int  kappa = kappa_from_lj(l, l+0.5);
  double sigma = w1*dirac_rr_nonneg(dirac_rr_xsec_dipole(eV, z, n, kappa), eV, n, kappa);
  if (l>0) {
    double w2 = 2.0*(l-0.5)+1.0;
    kappa = kappa_from_lj(l, l-0.5);
    sigma += w2*dirac_rr_nonneg(dirac_rr_xsec_dipole(eV, z, n, kappa), eV, n, kappa);
  }
  return sigma*eV;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief relativistic (Dirac) hydrogenic RR cross section in dipole approximation summed over subshells
 *
 * @param eV electron-ion collision energy in eV
 * @param z nuclear charge
 * @param fraction fractional populations for n,l subshells
 * @param use_fraction flags whether the above fractions should be used or not 
 * @param nmin principal quantum number of lowest shell
 * @param nmax prinicpal qauntum number of highest shell
 * @param lmin orbital angular momentum quantum number of lowest subshell
 *
 * @return RR cross section times electron-ion collision energy in cm2 eV
 */
double sigmarrdid(double eV, double z, vector<double> fraction, int use_fraction,
          int nmin, int nmax, int lmin)
{
  double sum = 0.0;
  for (int n=nmin; n<=nmax; n++) {
    int n12 = (n - 1) * n / 2;
    int l0 = (n==nmin) ? lmin : 0;
    for (int l=l0; l<n; l++) {
      double sigma = sigmarrdid(eV,z,n,l);
      if (use_fraction) sigma *=fraction[n12+l];
      sum += sigma;
    }
  }
  return sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief relativistic (Dirac) fully retarded hydrogenic RR cross section for nonrelativistic subshell n,l
 *
 * @param eV electron-ion collision energy in eV
 * @param z nuclear charge
 * @param n principal quantum number of initial subshell
 * @param l orbital angular momentum quantum number
 *
 * @return RR cross section times electron-ion collision energy in cm2 eV
 */
double sigmarrdir(double eV, double z, int n, int l)
{  
  double w1 = 2.0*(l+0.5)+1.0;
  int  kappa = kappa_from_lj(l, l+0.5);
  double sigma = w1*dirac_rr_nonneg(dirac_rr_xsec_retarded(eV, z, n, kappa), eV, n, kappa);
  if (l>0) {
    double w2 = 2.0*(l-0.5)+1.0;
    kappa = kappa_from_lj(l, l-0.5);
    sigma += w2*dirac_rr_nonneg(dirac_rr_xsec_retarded(eV, z, n, kappa), eV, n, kappa);
  }
  return sigma*eV;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief relativistic (Dirac) fully retarded hydrogenic RR cross section summed over subshells
 *
 * @param eV electron-ion collision energy in eV
 * @param z nuclear charge
 * @param fraction fractional populations for n,l subshells
 * @param use_fraction flags whether the above fractions should be used or not 
 * @param nmin principal quantum number of lowest shell
 * @param nmax prinicpal qauntum number of highest shell
 * @param lmin orbital angular momentum quantum number of lowest subshell
 *
 * @return RR cross section times electron-ion collision energy in cm2 eV
 */
double sigmarrdir(double eV, double z, vector<double> fraction, int use_fraction,
          int nmin, int nmax, int lmin)
{
  double sum = 0.0;
  for (int n=nmin; n<=nmax; n++) {
    int n12 = (n - 1) * n / 2;
    int l0 = (n==nmin) ? lmin : 0;
    for (int l=l0; l<n; l++) {
      double sigma = sigmarrdir(eV,z,n,l);
      if (use_fraction) sigma *=fraction[n12+l];
      sum += sigma;
    }
  }
  return sum;
}

//////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief low-energy Stobbe (Gaunt) correction factor to semiclassical RR cross section
 *
 * @param n principal qauntum number
 *
 * @param correction factor
 */
double stobbe(int n) {
  double lsum = 0.0;
  for (int l = 0; l < n; l++) {
    lsum += (2.0 * l + 1.0) * fosciBC(n, l);
  }
  return 3.0 * sqrt(3.0) * hydroconst::pi / 16 / n / n / n * lsum;
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief output of low-energy Stobbe correction factors to file
 */
void calc_stobbe(void) {
  int nmax;
  ofstream fout;

  cout << " Give nmax : ";
  cin >> nmax;
  fout.open("Stobbecorr.dat");
  for (int n = 1; n <= nmax; n++) // n starts at 1 (stobbe(0) divides by n^3=0)
    fout << n << " " << stobbe(n) << "\n";
  fout.close();
  cout << "\n Stobbe correction factors written to Stobbecorr.dat.\n\n";
}

//////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief precomputed low-energy Stobbe correction factors
 */
double kstobbe(int n) {
  const double stobbetab[100] = {
      0.7973014601263097, 0.876185137748039, 0.907508212449926,
      0.924743401639527,  0.93581638846388,  0.943608658751098,
      0.949431439764705,  0.953971669087505, 0.957626026642256,
      0.960640572468532,  0.963176553031456, 0.965344334771262,
      0.967222177180456,  0.96886722316249,  0.970322239758068,
      0.971619916920603,  0.972785701488747, 0.973839719665624,
      0.974798114035793,  0.975673993973899, 0.976478124451407,
      0.97721943395042,   0.97790539484623,  0.978542312298236,
      0.979135546461289,  0.979689685399302, 0.980208681071948,
      0.980695957327253,  0.981154496436549, 0.981586909013667,
      0.981995490945762,  0.982382270081938, 0.982749044779071,
      0.983097415924431,  0.983428813695211, 0.983744520043171,
      0.984045687685223,  0.984333356221232, 0.984608465876629,
      0.984871869270909,  0.985124341537149, 0.985366589057629,
      0.985599257032808,  0.985822936062623, 0.986038167888209,
      0.986245450417208,  0.986445242135483, 0.98663796599147,
      0.986824012825752,  0.98700374440719,  0.987177496127627,
      0.987345579399426,  0.987508283793644, 0.987665878951198,
      0.987818616294869,  0.987966730566065, 0.988110441207099,
      0.9882499536069,    0.988385460225757, 0.988517141612672,
      0.988645167327187,  0.988769696776052, 0.988890879973833,
      0.989008858235477,  0.989123764807859, 0.989235725446557,
      0.989344858943334,  0.989451277609236, 0.989555087717599,
      0.989656389910845,  0.989755279574496, 0.989851847181457,
      0.989946178609314,  0.990038355433108, 0.990128455195772,
      0.990216551658218,  0.990302715030843, 0.99038701218806,
      0.990469506867306,  0.990550259853813, 0.990629329152346,
      0.990706770146964,  0.990782635749768, 0.990856976539537,
      0.990929840891035,  0.99100127509572,  0.991071323474532,
      0.991140028483344,  0.991207430811654, 0.991273569474995,
      0.991338481901563,  0.991402204013441, 0.991464770302852,
      0.991526213903768,  0.991586566659212, 0.991645859184559,
      0.99170412092711,   0.991761380222187, 0.991817664345996,
      0.991872999565475};
  return (n >= 1 && n <= 100) ? stobbetab[n - 1] : 1.0;
}

//////////////////////////////////////////////////////////////////////
/**
 * @brief semi-classical RR cross section summed over subshells
 *
 * @param eV electron-ion collision energy in eV
 * @param z nuclear charge
 * @param fraction fractional populations for n,l subshells
 * @param use_fraction flags whether the above fractions should be used or not 
 * @param nmin principal quantum number of lowest shell
 * @param nmax prinicpal quantum number of highest shell
 * @param lmin orbital angular momentum quantum number of lowest subshell
 *
 * @return RR cross section times electron-ion collision energy in cm2 eV
 */
double sigmarrscl(double eV, double z, vector<double> fraction, int use_fraction,
           int nmin, int nmax, int lmin) {
  const double Pi = hydroconst::pi;
  const double Alpha = hydroconst::alpha;
  const double A0 = hydroconst::a0_cm;
  const double Rydberg = hydroconst::Ryd_eV;

  double z2 = z * z * Rydberg;
  double meanfraction = 1.0;
  int n12 = (nmin - 1) * nmin / 2;
  if (use_fraction) {
    meanfraction = 0.0;
    for (int l = 0; l < nmin; l++) {
      meanfraction += fraction[n12 + l] * (2 * l + 1);
    }
    meanfraction /= (nmin * nmin);
  } else {
    // calculate fraction of holes of nmin-shell from fraction of holes
    // of nmin,lmin-subshell stored in fraction[0]
    meanfraction =
        (nmin * nmin - lmin * lmin - (1 - fraction[0]) * (2 * lmin + 1)) /
        nmin / nmin;
  }
  double nsum = meanfraction * z2 / (1.0 + nmin * nmin * eV / z2) / nmin;
  if (z > 0)
    nsum *= kstobbe(nmin);
  for (int n = nmin + 1; n <= nmax; n++) {
    n12 = (n - 1) * n / 2;
    if (use_fraction) {
      meanfraction = 0.0;
      for (int l = 0; l < n; l++) {
        meanfraction += fraction[n12 + l] * (2 * l + 1);
      }
      meanfraction /= (n * n);
    } else {
      meanfraction = 1.0;
    }
    double nterm = meanfraction * z2 / (1.0 + n * n * eV / z2) / n;
    if (z > 0)
      nterm *= kstobbe(n);
    nsum += nterm;
  }
  return 32.0 * Pi * Alpha * Alpha * Alpha * A0 * A0 / 3.0 / sqrt(3.0) * nsum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief quantum mechanical (for n<=49) and semi-classical (for n>49) RR cross section summed over subshells
 *
 * @param eV electron-ion collision energy in eV
 * @param z nuclear charge
 * @param fraction fractional populations for n,l subshells
 * @param use_fraction flags whether the above fractions should be used or not 
 * @param nmin principal quantum number of lowest shell
 * @param nmax prinicpal quantum number of highest shell
 * @param lmin orbital angular momentum quantum number of lowest subshell
 *
 * @return RR cross section times electron-ion collision energy in cm2 eV
 */
double sigmarraqm(double eV, double z, vector<double> fraction,
                  int use_fraction, int nmin, int nmax, int lmin)
{
  int nqm = 49;
  nqm = nmax > nqm ? nqm : nmax;
  double sigma = sigmarrqm(eV, z, fraction, use_fraction, nmin, nqm, lmin);
  if (nmax > nqm) {
    double ratio =
        hydroconst::pi * 3.0 * sqrt(3.0) / 16.0 *
        sigmarrqm(eV, z, fraction, use_fraction, nqm + 1, nqm + 1, lmin) /
        stobbe(nqm + 1) /
        sigmarrscl(eV, z, fraction, use_fraction, nqm + 1, nqm + 1, lmin);
    sigma +=
        ratio * sigmarrscl(eV, z, fraction, use_fraction, nqm + 1, nmax, lmin);
  }
  return sigma;
}

//////////////////////////////////////////////////////////////////////
/**
 * @brief semi-classical (for n<n2) RR cross section summed over subshells plus integration over n>=n2
 *
 * for curiosity only, not used in production
 *
 * @param eV electron-ion collision energy in eV
 * @param z nuclear charge
 * @param fraction fractional populations for n,l subshells (not used)
 * @param n2 principal quantum number from where the integration starts
 * @param nmin principal quantum number of lowest shell
 * @param nmax prinicpal qautum number of highest shell
 * @param lmin orbital angluar momentum quantum number of lowest subshell
 *
 * @return RR cross section times electron-ion collision energy in cm2 eV
 */
double sigmarrscl2(double eV, double z, vector<double> fraction, int n2,
                   int nmin, int nmax, int lmin) { // e*sigma (eV cm^2)
  double sigma = sigmarrscl(eV, z, fraction, 0, nmin, n2 - 1, lmin);
  sigma +=
      2.1E-22 * hydroconst::Ryd_eV * z * z * log((double(nmax) / double(n2)));
  return sigma;
}

//////////////////////////////////////////////////////////////////////
/**
 * @brief plasma rate coefficient from semi classical cross section summed over subshells
 *
 * @param kT plasma temperature in eV
 * @param z nuclear charge
 * @param nmin minimum principal quantum number in sum
 * @param nmax maximum principal quantum number in sum
 * @param nele number of electrons in nmin-shell
 *
 * @return RR rate coefficient
 */
double alphaRRplasmaSCL(double kt, double z, int nmin, int nmax, int nele) {
  const double Pi = hydroconst::pi;
  const double Alpha = hydroconst::alpha;
  const double A0 = hydroconst::a0_cm;
  const double Rydberg = hydroconst::Ryd_eV;
  const double Clight = hydroconst::clight_cm_s;
  const double Melectron = hydroconst::mec2_eV;

  double z2r = Rydberg * z * z;
  double x = z2r / kt;
  double fac = 32.0 * Pi * Alpha * Alpha * Alpha * A0 * A0 / 3.0 / sqrt(3.0);
  fac *= 4.0 * Clight * z2r * z2r / sqrt(2 * Pi * kt * Melectron) / kt;
  double t = 1.0 - nele / (2.0 * nmin * nmin);
  double ratecoeff =
      fac / pow(double(nmin), 3.0) * e1exp(x / nmin / nmin) * t * kstobbe(nmin);
  for (int n = nmin + 1; n <= nmax; n++) {
    ratecoeff += kstobbe(n) * fac / pow(double(n), 3.0) * e1exp(x / n / n);
  }
  return ratecoeff;
}

//////////////////////////////////////////////////////////////////////
/**
 * @brief plasma cooling coefficient from semi classical cross section
 *
 * @param kT plasma temperature in eV
 * @param z nuclear charge
 * @param nmin minimum principal quantum number in sum
 * @param nmax maximum principal quantum number in sum
 * @param nele number of electrons in nmin-shell
 *
 * @return RR cooling coefficient
 */
double betaRRplasmaSCL(double kt, double z, int nmin, int nmax, int nele) {
  const double Pi = hydroconst::pi;
  const double Alpha = hydroconst::alpha;
  const double A0 = hydroconst::a0_cm;
  const double Rydberg = hydroconst::Ryd_eV;
  const double Clight = hydroconst::clight_cm_s;
  const double Melectron = hydroconst::mec2_eV;

  double z2r = Rydberg * z * z;
  double x = z2r / kt;
  double fac = 32.0 * Pi * Alpha * Alpha * Alpha * A0 * A0 / 3.0 / sqrt(3.0);
  fac *= 4.0 * Clight * z2r * z2r / sqrt(2 * Pi * kt * Melectron) / kt;
  double t = 1.0 - nele / (2.0 * nmin * nmin);
  double ratecoeff = fac / pow(double(nmin), 3.0) *
                     (1 - x / nmin / nmin * e1exp(x / nmin / nmin)) * t *
                     kstobbe(nmin);
  for (int n = nmin + 1; n <= nmax; n++) {
    ratecoeff += kstobbe(n) * fac / pow(double(n), 3.0) *
                 (1 - x / n / n * e1exp(x / n / n));
  }
  return ratecoeff;
}

//////////////////////////////////////////////////////////////////////
/**
 * @brief entry to interactive RR calculations
 */
void calc_sigmaRR(void) {
  double z, jmin, emin, emax, edelta, eV, t = 0.0, sigma = 0.0, IP = 0.0;
  int nint=1, l, nmin, lmin, nmax, nele=0, kappa = 0, choice = 0;
  int fdim = 2, use_fraction = 0;
  char answer;
  string filename, fnroot, header, line;
  ofstream fout, fout2;

  cout << "\n Select type of RR calculation ?\n";
  cout << "\n  1:  semiclassical calculation with Stobbe corrections (SCS)";
  cout << "\n  2:  semiclassical calculation without Stobbe corrections (SCL)";
  cout << "\n  3:  nonrelativistic quantum mechanical calculation, dipole approximation (NRD)";
  cout << "\n  4:  relativistic (Dirac) dipole approximation (DID)";
  cout << "\n  5:  relativistic (Dirac) fully retarded calculation (DIR)";
  cout << "\n  6:  NRD just one n,l for a range of energies";
  cout << "\n  7:  NRD, all n,l up to n=nmax at fixed energy";
  cout << "\n  8:  SCS up to nint-1, from nint to nmax integral dn (SCI)";
  cout << "\n  9:  photoionization from NRD RR (PI)";
  cout << "\n 10:  Ichihara & Eichler tabulated values (interpolated), one n,kappa (IE)";
  cout << "\n 11:  comparison: SCS vs. NRD";
  cout << "\n 12:  comparison: NRD vs. DID";
  cout << "\n 13:  comparison: DID vs. DIR";
  cout << "\n 14:  comparison: DIR vs. IE";
  while ((choice < 1) || (choice > 14)) {
    cout << "\n\n Make a choice ................................... : ";
    cin >> choice;
  }
    
  cout << "\n Give effective nuclear charge ....................: ";
  cin >> z;
  z = fabs(z);
  if (choice == 2) {
    z = -fabs(z); // flags that Stobbe corrections are not applied
  } 

  if (choice == 6) {
    cout << "\n Give n,l .........................................: ";
    cin >> nmin >> lmin;
  } else if (choice == 7) {
    cout << "\n Give maximum main quantum number nmax ............: ";
    cin >> nmax;
    cout << "\n Give electron energy in eV .......................: ";
    cin >> eV;
  } else if ((choice == 10) || (choice == 14)) {
    cout << "\n Give n, l, j .....................................: ";
    cin >> nmin >> lmin >> jmin;
    while (nmin>3) {
      cout << " ATTENTION: Ichihara & Eichler cross sections avalable only for n<=3!" << endl;
      cout << "\n Give n, l, j .....................................: ";
      cin >> nmin >> lmin >> jmin;
    }
    kappa = kappa_from_lj(lmin, jmin);
  } else if ((choice == 11) || (choice==12) || (choice==13)) {
    cout << "\n Give n, l ........................................: ";
    cin >> nmin >> lmin;
  } else {
    cout << "\n Give nmin, lmin, nmax ............................: ";
    cin >> nmin >> lmin >> nmax;
    cout << "\n Give number of electrons in minimum n,l shell ....: ";
    cin >> nele;
    t = 1.0 - nele / (4.0 * lmin + 2.0); /* fraction of holes in min n,l subshell */
    if (choice == 8) {
      use_fraction = 0;
      cout << "\n Give n where integration dn starts ...............: ";
      cin >> nint;
    } else if (choice == 9) {
      use_fraction = 0;
      cout << "\n Give first ionization potential in eV ............: ";
      cin >> IP;
    } else {
      cout << "\n Use surviving fractions from file ? (y/n) ........: ";
      cin >> answer;
      if ((answer == 'y') || (answer == 'Y')) {
        use_fraction = 1;
        fdim = nmax * (nmax + 1) / 2;
      }
    }
  }

  cout << "\n Give filename for output (*.sig) .................: ";
  cin >> fnroot;
  filename = fnroot + ".sig";
  fout.open(filename);

  if (choice == 7) {
    filename = fnroot + ".snl";
    fout2.open(filename);
    fout2 << "    n   lmax lmean lhalf interpolated\n";
    vector<double> es(nmax, 0.0);
    vector<double> esmean(nmax, 0.0);
    int choice2;
    cout << "\n Which kind of output ?";
    cout << "\n 1: sigma(n,l) in cm^2";
    cout << "\n 2: sigma(n,l) times energy in eV cm^2";
    cout << "\n 3: sigma(n,l)/sigma(n)";
    cout << "\n 4: sigma(n,l) relative to maximum value per n";
    cout << "\n Make a choice ...................................: ";
    cin >> choice2;
    for (int n = 1; n <= nmax; n++) {
      double esmax = 0.0, estot = 0.0, esout;
      int lmax = 0;
      for (l = 0; l < n; l++) {
        es[l] = sigmarrqm(eV, z, n, l);
        estot += es[l];
        esmean[l] = estot / (l + 1.0);
        if (es[l] > esmax) {
          lmax = l;
          esmax = es[l];
        }
      }
      switch (choice2) {
      case 2:
        for (l = 0; l < n; l++) {
          esout = es[l] > 1e-99 ? es[l] : 1e-99;
          fout << setw(12) << setprecision(6) << esout;
        }
        break;
      case 3:
        for (l = 0; l < n; l++) {
          fout << setw(12) << setprecision(6) << es[l] / estot;
        }
        break;
      case 4:
        for (l = 0; l < n; l++) {
          fout << setw(12) << setprecision(6) << es[l] / esmax;
        }
        break;
      default:
        for (l = 0; l < n; l++) {
          esout = es[l] / eV > 1e-99 ? es[l] / eV : 1e-99;
          fout << setw(12) << setprecision(6) << esout;
        }
        break;
      }
      for (int k = n; k < nmax; k++) {
        fout << setw(12) << setprecision(6) << 1e-99;
      }
      fout << "\n";

      int lhalf = 0, lmean = 0;
      for (int ll = n - 1; ll >= 0; ll--) {
        if ((lhalf == 0) && (es[ll] > 0.5 * es[lmax])) {
          lhalf = ll;
        }
        if ((lmean == 0) && (es[ll] > 0.5 * esmean[ll])) {
          lmean = ll;
        }
      }

      double ilhalf = lhalf > 0 ? lhalf + (0.5 * es[lmax] - es[lhalf]) /
                                              (es[lhalf + 1] - es[lhalf])
                                : 0;
      fout2 << setw(5) << n << " " << setw(5) << lmax << " " << setw(5) << lmean
            << " " << setw(5) << lhalf << " " << setw(12) << setprecision(5)
            << ilhalf << "\n";
    } // end for(n...)
    fout.close();
    fout2.close();
    cout << "\n Quantum mechanical RR (NRD) cross sections at E = " << setw(10)
         << setprecision(3) << eV << " eV";
    cout << "\n written in matrix form to file " << fnroot << ".sig.\n";
    return;
  }

  vector<double> fraction(fdim, 1.0);
  if (use_fraction) {
    for (int k = 0; k < fdim; k++) {
      fraction[k] = 1.0;
    }
    for (l = 0; l < lmin; l++) {
      fraction[(nmin - 1) * nmin / 2 + l] = 0.0;
    }
    fraction[(nmin - 1) * nmin / 2 + lmin] = t;
    nmax = readfraction(fraction, header, nmax);
    cout << "\n surviving fractions " << header << "\n";
    fout << "surviving fractions " << header << "\n";
  } else {
    fraction[0] = t;
  }

  switch (choice) {
  case 1:
    fout << "semiclassical RR cross section";
    fout << "  with Stobbe correction\n";
    fout << "SCS: z=" << fixed << setprecision(2) << z << ", nmin=" << setw(3)
         << nmin << ", lmin=" << setw(3) << lmin;
    fout << " nele=" << setw(3) << nele << ", nmax=" << setw(3) << nmax << "\n";
    break;
  case 2:
    fout << "semiclassical RR cross section";
    fout << "  without Stobbe correction\n";
    fout << "SCL: z=" << fixed << setprecision(2) << z << ", nmin=" << setw(3)
         << nmin << ", lmin=" << setw(3) << lmin;
    fout << " nele=" << setw(3) << nele << ", nmax=" << setw(3) << nmax << "\n";
    break;
  case 3:
    fout << "nonrelativistic quantum mechanical RR cross section";
    fout << " (dipole approximation)\n";
    fout << "NRD: z=" << fixed << setprecision(2) << z << ", nmin=" << setw(3)
         << nmin << ", lmin=" << setw(3) << lmin;
    fout << " nele=" << setw(3) << nele << ", nmax=" << setw(3) << nmax << "\n";
    break;
  case 4:
    fout << " relativistic (Dirac) RR cross section (dipole approximation)";
    fout << " \n";
    fout << "DID: z=" << fixed << setprecision(2) << z << ", n=" << setw(3)
         << nmin << ", l=" << setw(3) << lmin << ", j=" << setw(3) << jmin << "\n";
    break;
  case 5:
    fout << " relativistic fully retarded RR cross section";
    fout << " \n";
    fout << "DIR: z=" << fixed << setprecision(2) << z << ", n=" << setw(3)
         << nmin << ", l=" << setw(3) << lmin << ", j=" << setw(3) << jmin << "\n";
    break;
  case 6:
    fout << "nonrelativistic quantum mechanical RR cross section";
    fout << " (dipole approximation)\n";
    fout << "NRD: z=" << fixed << setprecision(2) << z << ", n=" << setw(3)
         << nmin << ", l=" << setw(3) << lmin << "\n";
    break;
  case 8:
    fout << "semiclassical RR cross section";
    fout << "  with high-n integration\n";
    fout << "SCI: z=" << fixed << setprecision(2) << z << ", n=" << setw(3)
         << nmin << ", l=" << setw(3) << lmin << "\n";
    fout << " nele=" << setw(3) << nele << ", nint=" << setw(3) << nint
         << ", nmax=" << setw(3) << nmax << "\n";
    break;
  case 9:
    fout << "non relativistic PI cross section";
    fout << "  (dipole approximation) \n";
    fout << "PI: z=" << fixed << setprecision(2) << z << ", n=" << setw(3)
         << nmin << ", l=" << setw(3) << lmin << "\n";
    fout << " nele=" << setw(3) << nele << ", nmax=" << setw(3) << nmax
         << ", IP=" << fixed << setprecision(3) << IP << "\n";
    break;
  case 10:
    fout << "Ichihara & Eichler tabulated RR cross section (interpolated)\n";
    fout << "IE: z=" << fixed << setprecision(2) << z << ", n=" << setw(3)
         << nmin << ", l=" << setw(3) << lmin << ", j=" << setw(3) << jmin << "\n";
    break;
  case 11:
    fout << "comparison: semiclassical (SCL) vs. nonrelativistic dipole (NRD)\n";
    fout << "SCL/NRD: z=" << fixed << setprecision(2) << z << ", n=" << setw(3)
         << nmin << ", l=" << setw(3) << lmin << "\n";
    break;
  case 12:
    fout << "comparison: nonrelativistic dipole (NRD) vs. relativistic dipole (DID)\n";
    fout << "NRD/DID: z=" << fixed << setprecision(2) << z << ", n=" << setw(3)
         << nmin << ", l=" << setw(3) << lmin << "\n";
    break;
  case 13:
    fout << "comparison: Dirac dipole (DID) vs. Dirac fully retarded (DIR) \n";
    fout << "DID/DIR: z=" << fixed << setprecision(2) << z << ", n=" << setw(3)
         << nmin << ", l=" << setw(3) << lmin << "\n";
    break;
  case 14:
    fout << "comparison: Dirac fully retarded (DIR) vs. Ichihara & Eichler (IE)\n";
    fout << "DIR/IE: z=" << fixed << setprecision(2) << z << ", n=" << setw(3)
         << nmin << ", l=" << setw(3) << lmin << ", j=" << setw(3) << jmin << "\n";
    break;
  }

  int steps;
  ifstream fenergy;
  line = "0";
  cout << "\n Read energies from file? Give filename (0=no file): ";
  cin >> filename;
  int read_energy = (filename != "0");

  if (read_energy) {
    fenergy.open(filename);
    if (!fenergy.is_open()) {
      read_energy = 0;
      cout << "\n file " << filename << " not found!";
    } else { // count the number of energies given in the data file
      steps = 0;
      while (getline(fenergy, line))  steps++;
      steps--;
      fenergy.close();
    }
  }
  if (!read_energy) {
    cout << "\n Give energy range (emin,emax,delta) ..............: ";
    cin >> emin >> emax >> edelta;
    steps = int((emax - emin) / edelta);
  }

  if (choice == 11) {
    cout << "\n     Ecm (eV)     sSCL (cm^2)      sNRD (cm^2)     sSCL/sNRD\n";
    fout << "     Ecm (eV)     sSCL (cm^2)      sNRD (cm^2)     sSCL/sNRD\n";
  } else if (choice == 12) {
    cout << "\n     Ecm (eV)     sNRD (cm^2)      sDID (cm^2)     sNRD/sDID\n";
    fout << "     Ecm (eV)     sNRD (cm^2)      sDID (cm^2)     sNRD/sDID\n";
  } else if (choice == 13) {
    cout << "\n     Ecm (eV)     sDID (cm^2)      sDIR (cm^2)     sDID/sDIR\n";
    fout << "     Ecm (eV)     sDID (cm^2)      sDIR (cm^2)     sDID/sDIR\n";
  } else if (choice == 14) {
    cout << "\n     Ecm (eV)     sDIR (cm^2)      sIE (cm^2)      sDIR/sIE\n";
    fout << "     Ecm (eV)     sDIR (cm^2)      sIE (cm^2)      sDIR/sIE\n";
  } else {
    cout << "\n     Ecm (eV)        sigma (cm^2)    E*sigma (eV cm^2)\n";
    fout << "     Ecm             sigma           E*sigma\n";
    fout << "     (eV)            (cm^2)          (eV cm^2)\n";
  }

  if (read_energy)
    fenergy.open(filename);

  for (int i = 0; i <= steps; i++) {
    if (read_energy) {
      getline(fenergy, line);
      istringstream ss(line);
      ss >> eV;
    } else {
      eV = emin + i * edelta;
      if (emin < 0)
        eV = exp(log(10.0) * eV);
    }

    switch (choice) {
    case 1:
      sigma = sigmarrscl(eV, z, fraction, use_fraction, nmin, nmax, lmin);
      break;
    case 2:
      sigma = sigmarrscl(eV, z, fraction, use_fraction, nmin, nmax, lmin);
      break;
    case 3:
      sigma = sigmarrqm(eV, z, fraction, use_fraction, nmin, nmax, lmin);
      break;
    case 4:
      sigma = sigmarrdid(eV, z, fraction, use_fraction, nmin, nmax, lmin);
      break;
    case 5:
      sigma = sigmarrdir(eV, z, fraction, use_fraction, nmin, nmax, lmin);
      break;
    case 8:
      sigma = sigmarrscl2(eV, z, fraction, nint, nmin, nmax, lmin);
      break;
    case 9:
      sigma = sigmaPIqm(eV, z, fraction, use_fraction, nmin, nmax, lmin, IP);
      break;
    case 10: {
      sigma = sigma_IchiharaEichlerRR(eV, (int)z, nmin, lmin, jmin) * eV;
      break;
    }
    case 11: {
      double sigma_qmd = 0.0, sigma_scl = 0.0;
      if (eV>0.0) {
	 sigma_qmd = sigmarrqm(eV, z, fraction, use_fraction, nmin, nmin, lmin)/eV;
	 sigma_scl = sigmarrscl(eV, z, fraction, use_fraction, nmin, nmin, lmin)/eV;
      }
      double ratio = sigma_qmd > 0.0 ? sigma_scl / sigma_qmd : 0.0;
      cout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
	   << eV << " " << setw(15) << uppercase << defaultfloat
	   << setprecision(5) << sigma_scl << " " << setw(15) << uppercase
	   << defaultfloat << setprecision(5) << sigma_qmd << " " << setw(12)
	   << uppercase << defaultfloat << setprecision(5) << ratio << "\n";
      fout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
         << eV << " " << setw(15) << uppercase << defaultfloat
	   << setprecision(5) << sigma_scl << " " << setw(15) << uppercase
         << defaultfloat << setprecision(5) << sigma_qmd << " " << setw(12)
         << uppercase << defaultfloat << setprecision(5) << ratio << "\n";
      break;
    }
    case 12: {
      double sigma_qmd = 0.0, sigma_did= 0.0;
      if (eV>0.0) {
	sigma_qmd = sigmarrqm(eV, z, nmin, lmin)/eV;
	sigma_did = sigmarrdid(eV, z, nmin, lmin)/eV;
      }
      double ratio = sigma_qmd > 0.0 ? sigma_qmd / sigma_did : 0.0;
      cout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
	   << eV << " " << setw(15) << uppercase << defaultfloat
	   << setprecision(5) << sigma_qmd << " " << setw(15) << uppercase
	   << defaultfloat << setprecision(5) << sigma_did << " " << setw(12)
	   << uppercase << defaultfloat << setprecision(5) << ratio  << "\n";
      fout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
         << eV << " " << setw(15) << uppercase << defaultfloat
	   << setprecision(5) << sigma_qmd << " " << setw(15) << uppercase
         << defaultfloat << setprecision(5) << sigma_did << " " << setw(12)
         << uppercase << defaultfloat << setprecision(5) << ratio << "\n";
      break;
    }
    case 13: {
      double sigma_did= 0.0, sigma_dir=0.0;
      if (eV>0.0) {
	sigma_did = sigmarrdid(eV, z, nmin, lmin)/eV;
	sigma_dir = sigmarrdir(eV, z, nmin, lmin)/eV;
      }
      double ratio = sigma_dir > 0.0 ? sigma_did / sigma_dir : 0.0;
      cout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
	   << eV << " " << setw(15) << uppercase << defaultfloat
	   << setprecision(5) << sigma_did << " " << setw(15) << uppercase
	   << defaultfloat << setprecision(5) << sigma_dir << " " << setw(12)
	   << uppercase << defaultfloat << setprecision(5) << ratio  << "\n";
      fout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
         << eV << " " << setw(15) << uppercase << defaultfloat
	   << setprecision(5) << sigma_did << " " << setw(15) << uppercase
         << defaultfloat << setprecision(5) << sigma_dir << " " << setw(12)
         << uppercase << defaultfloat << setprecision(5) << ratio << "\n";
      break;
    }
    case 14: {
      double sigma_dir = dirac_rr_nonneg(dirac_rr_xsec_retarded(eV, z, nmin, kappa), eV, nmin, kappa);
      double sigma_ie = sigma_IchiharaEichlerRR(eV, (int)z, nmin, lmin, jmin);
      double ratio = sigma_ie > 0.0 ? sigma_dir / sigma_ie : 0.0;
      cout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
           << eV << " " << setw(15) << uppercase << defaultfloat
           << setprecision(5) << sigma_dir << " " << setw(15) << uppercase
           << defaultfloat << setprecision(5) << sigma_ie << " " << setw(12)
           << uppercase << defaultfloat << setprecision(5) << ratio << "\n";
      fout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
           << eV << " " << setw(15) << uppercase << defaultfloat
           << setprecision(5) << sigma_dir << " " << setw(15) << uppercase
           << defaultfloat << setprecision(5) << sigma_ie << " " << setw(12)
           << uppercase << defaultfloat << setprecision(5) << ratio << "\n";
      break;
    }
    default:
      sigma = sigmarrqm(eV, z, nmin, lmin);
      break;
    }
    if (choice < 10) {
      cout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
	   << eV << "     " << setw(15) << uppercase << defaultfloat
	   << setprecision(5) << sigma / eV << "      " << setw(15) << uppercase
	   << defaultfloat << setprecision(5) << sigma << "\n";
      fout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
	   << eV << "     " << setw(15) << uppercase << defaultfloat
	   << setprecision(5) << sigma / eV << "      " << setw(15) << uppercase
	   << defaultfloat << setprecision(5) << sigma << "\n";
    }
  }
  if (read_energy)
    fenergy.close();
  fout.close();
}


//////////////////////////////////////////////////////////////////////

void calcAlphaRRplasma(void) {
  double z;
  int nmin, nmax, nele;
  string filename, line;
  ifstream ftemp;
  ofstream fout;

  cout << "\n Give z, nmin, nmax ........................................: ";
  cin >> z >> nmin >> nmax;
  cout << "\n Give number of electrons already in nmin-shell ............: ";
  cin >> nele;

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

  double kt, ktmin, ktmax, ktdelta;
  if (read_temp) {
    ftemp.open(filename);
  } else {
    cout << "\n Give log min, max, delta kT in eV  ........................: ";
    cin >> ktmin >> ktmax >> ktdelta;
    steps = int((ktmax - ktmin) / ktdelta);
  }

  cout << "\n Give filename for output ..................................: ";
  cin >> filename;
  fout.open(filename);
  fout << "z=" << fixed << setprecision(1) << z << ", nmin=" << setw(2) << nmin
       << ", nele=" << setw(2) << nele << ", nmax=" << setw(3) << nmax << "\n";

  cout << "\n      kT (eV)      alpha (cm^3/s)        kT^3/2*alpha         "
          "beta (cm^3/s)      kT^(5/2)*beta\n";
  fout << "      kT (eV)      alpha (cm^3/s)        kT^3/2*alpha          beta "
          "(cm^3/s)     kT^5/2*beta\n";

  for (int i = 0; i <= steps; i++) {
    if (read_temp) {
      getline(ftemp, line);
      istringstream ss(line);
      ss >> kt;
    } else {
      kt = ktmin + i * ktdelta;
      kt = exp(log(10.0) * kt);
    }
    double alpha = alphaRRplasmaSCL(kt, z, nmin, nmax, nele);
    double beta = betaRRplasmaSCL(kt, z, nmin, nmax, nele);
    double akt = alpha * kt * sqrt(kt);
    double bkt = beta * kt * kt * sqrt(kt);
    cout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
         << kt << "     " << setw(15) << uppercase << defaultfloat
         << setprecision(5) << alpha << "      " << setw(15) << uppercase
         << defaultfloat << setprecision(5) << akt << "     " << setw(15)
         << uppercase << defaultfloat << setprecision(5) << beta << "     "
         << setw(15) << uppercase << defaultfloat << setprecision(5) << bkt
         << "\n";
    fout << " " << setw(12) << uppercase << defaultfloat << setprecision(5)
         << kt << "     " << setw(15) << uppercase << defaultfloat
         << setprecision(5) << alpha << "      " << setw(15) << uppercase
         << defaultfloat << setprecision(5) << akt << "     " << setw(15)
         << uppercase << defaultfloat << setprecision(5) << beta << "     "
         << setw(15) << uppercase << defaultfloat << setprecision(5) << bkt
         << "\n";
  }
  fout.close();
}

//////////////////////////////////////////////////////////////////////
