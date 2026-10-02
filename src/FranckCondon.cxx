/**
 * @file FranckCondon.cxx
 *
 * @brief Calculation of Franck-Condon factors of diatomic molecules
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: FranckCondon.cxx 614 2019-03-29 14:01:00Z iamp $
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */
#include "gaussint.h"
#include "hydroconst.h"
#include "matrix.h"
#include "peakfunctions.h"
#include "buildinfo.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>

#undef useBOOST
#if __has_include(<boost/multiprecision/cpp_dec_float.hpp>)
#if __has_include(<boost/math/special_functions/gamma.hpp>)
#if __has_include(<boost/math/quadrature/gauss_kronrod.hpp>)
#if __has_include(<boost/math/quadrature/exp_sinh.hpp>)
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::math;
using namespace boost::multiprecision;
typedef cpp_dec_float_100 mp_double;
#define useBOOST
#endif
#endif
#endif
#endif

#ifndef useBOOST
typedef long double mp_double;
#endif

using namespace std;
using hydroconst::pi;

// the following switch can be activated for testing the numerical integration
// #define testMorseJ

// coefficients for gaussian quadrature rules
const int FClegendre_n = 96, FClaguerre_n = 15;
vector<mp_double> FClegendre_x, FClegendre_w;
vector<mp_double> FClaguerre_x, FClaguerre_w;

void FCsetup_integration(void) {
  vector<double> legendre_x, legendre_w;
  gauss_coef(FClegendre_n, legendre_x, legendre_w); // see gaussint.cxx
  FClegendre_x.resize(FClegendre_n);
  FClegendre_w.resize(FClegendre_n);
  for (int n = 0; n < FClegendre_n; n++) {
    FClegendre_x[n] = (mp_double)legendre_x[n];
    FClegendre_w[n] = (mp_double)legendre_w[n];
  }

  vector<double> laguerre_x, laguerre_w;
  laguer_coef(FClaguerre_n, laguerre_x, laguerre_w); // see gaussint.cxx
  FClaguerre_x.resize(FClaguerre_n);
  FClaguerre_w.resize(FClaguerre_n);
  for (int n = 0; n < FClaguerre_n; n++) {
    FClaguerre_x[n] = (mp_double)laguerre_x[n];
    FClaguerre_w[n] = (mp_double)laguerre_w[n];
  }
}

//**************************************************************************************************************
#ifdef useBOOST
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Integral appearing the in the calculation of the Franc-Condon factors
 * for two displaced Morse potentials
 *
 * using adaptive quadrature from the BOOST library
 *
 * Eq. 22 of Lopez et al., International Journal of Quantum Chemistry, Vol 88,
 * 280295 (2002)
 */
mp_double boost_MorseFCintegralJ(mp_double a, mp_double k1, mp_double k2,
                                 mp_double y1, mp_double y2) {
  mp_double error;
  mp_double kak = k1 + a * k2 - 1.0;

  // the integrand is defined as an anonymous function (also called lambda
  // function) see
  // https://en.wikipedia.org/wiki/Anonymous_function#C++_(since_C++11)
  auto integrand = [kak, y1, y2, a](mp_double u) -> mp_double {
    return exp(log(u) * kak - 0.5 * (y1 * u + y2 * pow(u, a)));
  };

  mp_double result =
      boost::math::quadrature::gauss_kronrod<mp_double, 61>::integrate(
          integrand, (mp_double)0.0, numeric_limits<mp_double>::infinity(), 5,
          1e-24, &error);

  // cout << scientific << result << "   " << error << "   " << error*L1 <<
  // endl;
  return result;
}
#endif
//**************************************************************************************************************

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Integral appearing the in the calculation of the Franc-Condon factors
 * for two displaced Morse potentials
 *
 * numerical integration by gaussian quadrature
 *
 * Eq. 22 of Lopez et al., International Journal of Quantum Chemistry, Vol 88,
 * 280295 (2002)
 */
mp_double MorseFCintegralJ(mp_double a, mp_double k1, mp_double k2,
                           mp_double y1, mp_double y2) {
#ifdef testMorseJ
  a = 1.0;
#endif
  mp_double kak = k1 + a * k2 - 1.0;

  // search for maximum of integrand, i.e, the zero of its derivative
  mp_double ulo = 0.0;            // initial lower bound of the search interval
  mp_double uhi = 2.0 * kak / y1; // initial upper bound of the search interval
  mp_double umax = 0.0;
  int count = 0;
  while ((fabs(uhi - ulo) > 1.0E-3) && (count < 100)) {
    count++; // typically eleven iterations were required in the test cases
    umax = 0.5 * (ulo + uhi); // bisection of the search interval
    mp_double test =
        kak - 0.5 * y1 * umax -
        0.5 * a * y2 *
            pow(umax,
                a); // decisive factor from the derivative of the integrand
    if (test < 0.0)
      uhi = umax;
    else
      ulo = umax;
  }
  // now 'umax' is the position of the maximum of the integrand

  // subdivision of integration intervals chosen such that test cases are
  // reproduced as close as achievable
  const int intervals = 4;
  mp_double u[intervals + 1];
  u[0] = 0.0;
  u[1] = 0.7 * umax;
  u[2] = 1.3 * umax;
  u[3] = 7.0 * umax;
  u[4] = u[3] + 100.0 / (y1 + y2); // cutoff

  mp_double scalefactor =
      exp(kak * log(umax) -
          0.5 * (y1 * umax + y2 * pow(umax, a))); // max value of the integrand

  // Gauss-Legendre integration on each finite subinterval
  mp_double sum = 0.0, result = 0.0;
  for (int i = 1; i <= intervals; i++) {
    sum = 0.0;
    for (int n = FClegendre_n - 1; n >= 0; n--) {
      mp_double uu =
          0.5 * ((u[i] - u[i - 1]) * FClegendre_x[n] + u[i] + u[i - 1]);
      mp_double y = FClegendre_w[n] *
                    exp(kak * log(uu) - 0.5 * (y1 * uu + y2 * pow(uu, a))) /
                    scalefactor;
      sum += y;
    }
    result += sum * 0.5 * (u[i] - u[i - 1]);
  }

  // Gauss-Laguerre integration from cutoff to infinity
  sum = 0.0;
  for (int n = 0; n < FClaguerre_n; n++) {
    mp_double uu = FClaguerre_x[n] + u[intervals];
    mp_double y = FClaguerre_w[n] *
                  exp(kak * log(uu) + FClaguerre_x[n] -
                      0.5 * (y1 * uu + y2 * pow(uu, a))) /
                  scalefactor;
    sum += y;
  }
  result += sum;
  result *= scalefactor;
#ifdef testMorseJ
  mp_double result1 =
      exp(lgamma(k1 + k2) - (k1 + k2) * log(0.5 * (y1 + y2))); // Eq. 25 for a=1
  cout << kak << "  " << result << " " << result / result1 - 1.0 << endl;
#endif
  return result;
}

///////////////////////////////////////////////////////////////////////////////
/**
 * @brief Franck-Condon factor for two displaced Morse-Potentials
 *
 * Lopez et al., International Journal of Quantum Chemistry, Vol 88, 280295
 * (2002).
 */
double MorseFC(int n1, int n2, double omega1, double omega2, double omegachi1,
               double omegachi2, double reduced_mass, double x0,
               int test_mode = 0) {

  // Eq. 20 of Lopez et al., International Journal of Quantum Chemistry, Vol 88,
  // 280295 (2002)

  if (test_mode > 0)
    omegachi2 = omegachi1;                         // leads to a=1 below
  const mp_double hbarc = hydroconst::hbarc_eV_nm; // hbar*c in eV nm
  mp_double j1 = 0.5 * (omega1 / omegachi1 - 1.0);
  mp_double j2 = 0.5 * (omega2 / omegachi2 - 1.0);
  mp_double beta1 = sqrt(2.0 * reduced_mass * omegachi1) / hbarc; // in 1/nm
  mp_double beta2 = sqrt(2.0 * reduced_mass * omegachi2) / hbarc; // in 1/nm

  mp_double a = beta2 / beta1;
  mp_double y1 = (2.0 * j1 + 1.0) * exp(-0.5 * beta1 * x0);
  mp_double y2 = (2.0 * j2 + 1.0) * exp(0.5 * beta2 * x0);
  mp_double j1n1 = j1 - n1;
  mp_double j2n2 = j2 - n2;
  mp_double twoj1n1 = j1 + j1n1 + 1.0;
  mp_double twoj2n2 = j2 + j2n2 + 1.0;
  mp_double lgamma1 = lgamma(twoj1n1);
  mp_double lgamma2 = lgamma(twoj2n2);

  mp_double term0 =
      log(2.0) + 0.5 * (log(a) + lgamma(n1 + 1.0) + lgamma(n2 + 1.0) - lgamma1 -
                        lgamma2 + log(j1n1) + log(j2n2));
  mp_double sum = 0.0, c = 0.0;
  for (int l1 = 0; l1 <= n1; l1++) {
    for (int l2 = 0; l2 <= n2; l2++) {
      int n1l1 = n1 - l1;
      int n2l2 = n2 - l2;
      mp_double binomial1 =
          lgamma1 - lgamma(twoj1n1 - n1l1) - lgamma(n1l1 + 1.0);
      mp_double binomial2 =
          lgamma2 - lgamma(twoj2n2 - n2l2) - lgamma(n2l2 + 1.0);
      mp_double term1 = (j1n1 + l1) * log(y1) + (j2n2 + l2) * log(y2) -
                        lgamma(l1 + 1.0) - lgamma(l2 + 1.0) + binomial1 +
                        binomial2;
      // cout << term << endl;
      int sign = pow(-1, l1 + l2);
      mp_double k1 = j1n1 + l1;
      mp_double k2 = j2n2 + l2;
      mp_double lnJ;
      if (test_mode == 2) {
        lnJ = lgamma(k1 + k2) -
              (k1 + k2) * log(0.5 * (y1 + y2)); // Eq. 25 for a=1
      } else {
#ifdef useBOOST
        lnJ = log(boost_MorseFCintegralJ(a, k1, k2, y1, y2));
#else
        lnJ = log(MorseFCintegralJ(a, k1, k2, y1, y2));
#endif
      }
      mp_double term = sign * exp(term0 + term1 + lnJ);
      // cout << term << endl;
      mp_double t =
          sum + term; // summation according to the improved Kahan-Babuška
                      // algorithm to keep track of round-off errors
      if (fabs(sum) > fabs(term)) {
        c += (sum - t) + term;
      } else {
        c += (term - t) + sum;
      }
      sum = t;
    }
  }
  // cout << scientific << c << endl;
  return (double)(sum + c);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculation of Franck-Condon factors for harmonic and Morse potentials
 */
int FranckCondon(void) {

  // harmonic oscillator frequencies in eV (hbar*omega)
  double omega1 = 0.5, omega2 = 0.5, omegachi1 = 0.0, omegachi2 = 0.0;

  double Re1; // equilibrium distance of lower oscillator
  double Re2; // equilibrium distance of upper oscillator
  // spacing between harmonic oscillators in nm
  double x0 =
      0.093; // = Re2-Re1 as in Lopez et al., Frank et al. use the opposite sign

  // atomic mass numbers of atomic constituents of the diatomic molecule
  double A1 = 6, A2 = 6, Areduced;

  // max vibrational quantum numbers;
  int v1max = 4, v2max = 4;

  char answer;
  int selected_test = 0;

  cout << " Calculation of Franck condon factors for transitions between two "
          "potentials of a diatomic molecule.";
  cout << endl
       << " Harmonic and  Morse potentials are assumed." << endl
       << endl;
  cout << endl << " Arbitrary input values or test cases ? (a/t) ..........: ";
  cin >> answer;
  if ((answer == 't') || (answer == 'T')) {
    A1 = 1.0;
    A2 = 1.0; // reduced mass is 1/2!
    cout << endl
         << " The calcuation of the Franck Condon factors for the Morse "
            "potential follow";
    cout << endl
         << " the prescription of Lopez et al, Int. J. Quant. Chem. 88 (2002) "
            "280295,";
    cout << endl
         << " which requires a numerical integration. The following test cases "
            "can be";
    cout << endl
         << " used to check the accuracy of the numerical integration: ";
    cout << endl << "   1) Tab. IV of Lopez et al.";
    cout << endl << "   2) Tab. V of Lopez et al.";
    cout << endl << "   3) Tab. VI of Lopez et al.";
    cout << endl << "   4) Tab. VII of Lopez et al.";
    cout << endl
         << " Which test case to use? (1/2/3/4) .....................: ";
    cin >> selected_test;
  }
  if (selected_test ==
      1) { // should reproduce Tab. IV of Lopez et al., International Journal of
           // Quantum Chemistry, Vol 88, 280295 (2002)
    v1max = 6;
    v2max = 6;
    x0 = -(5.0535 - 5.876);
    double j2 = 69.494, j1 = 48.106;
    double beta2 = 0.326, beta1 = 0.459;
    cout << endl
         << " Test case IV: j1 = " << j1 << ", j2 = " << j2
         << ", beta1 = " << beta1 << "/a0, beta2 = " << beta2
         << "/a0, R2-R1 = " << x0 << " a0";
    omegachi1 = pow(hydroconst::hbarc_eV_nm * beta1 / hydroconst::a0_nm, 2) /
                hydroconst::muc2_eV;
    omegachi2 = pow(hydroconst::hbarc_eV_nm * beta2 / hydroconst::a0_nm, 2) /
                hydroconst::muc2_eV;
    omega1 = (2.0 * j1 + 1.0) * omegachi1;
    omega2 = (2.0 * j2 + 1.0) * omegachi2;
    Areduced = 0.5 *
               pow(beta1 * hydroconst::hbarc_eV_nm / hydroconst::a0_nm, 2) /
               omegachi1 / hydroconst::muc2_eV;
    x0 *= hydroconst::a0_nm;
  } else if (selected_test ==
             2) { // should reproduce Tab. V of Lopez et al., International
                  // Journal of Quantum Chemistry, Vol 88, 280295 (2002)
    v1max = 5;
    v2max = 6;
    x0 = -(5.0535 - 5.550);
    double j2 = 62.05, j1 = 48.104;
    double beta2 = 0.354, beta1 = 0.459;
    cout << endl
         << " Test case V: j1 = " << j1 << ", j2 = " << j2
         << ", beta1 = " << beta1 << "/a0, beta2 = " << beta2
         << "/a0, R2-R1 = " << x0 << " a0";
    omegachi1 = pow(hydroconst::hbarc_eV_nm * beta1 / hydroconst::a0_nm, 2) /
                hydroconst::muc2_eV;
    omegachi2 = pow(hydroconst::hbarc_eV_nm * beta2 / hydroconst::a0_nm, 2) /
                hydroconst::muc2_eV;
    omega1 = (2.0 * j1 + 1.0) * omegachi1;
    omega2 = (2.0 * j2 + 1.0) * omegachi2;
    Areduced = 0.5 *
               pow(beta1 * hydroconst::hbarc_eV_nm / hydroconst::a0_nm, 2) /
               omegachi1 / hydroconst::muc2_eV;
    x0 *= hydroconst::a0_nm;
  } else if (selected_test ==
             3) { // should reproduce Tab. VI of Lopez et al., International
                  // Journal of Quantum Chemistry, Vol 88, 280295 (2002)
    v1max = 3;
    v2max = 3;
    x0 = -(2.2152 - 2.3308);
    double j2 = 60.4084, j1 = 73.5602;
    double beta2 = 1.2167, beta1 = 1.227;
    cout << endl
         << " Test case VI: j1 = " << j1 << ", j2 = " << j2
         << ", beta1 = " << beta1 << "/a0, beta2 = " << beta2
         << "/a0, R2-R1 = " << x0 << " a0";
    omegachi1 = pow(hydroconst::hbarc_eV_nm * beta1 / hydroconst::a0_nm, 2) /
                hydroconst::muc2_eV;
    omegachi2 = pow(hydroconst::hbarc_eV_nm * beta2 / hydroconst::a0_nm, 2) /
                hydroconst::muc2_eV;
    omega1 = (2.0 * j1 + 1.0) * omegachi1;
    omega2 = (2.0 * j2 + 1.0) * omegachi2;
    Areduced = 0.5 *
               pow(beta1 * hydroconst::hbarc_eV_nm / hydroconst::a0_nm, 2) /
               omegachi1 / hydroconst::muc2_eV;
    x0 *= hydroconst::a0_nm;
  } else if (selected_test ==
             4) { // should reproduce Tab. VII of Lopez et al., International
                  // Journal of Quantum Chemistry, Vol 88, 280295 (2002)
    v1max = 6;
    v2max = 6;
    x0 = -(2.420 - 2.287);
    double j2 = 58.55, j1 = 50.994;
    double beta2 = 1.301, beta1 = 1.280;
    cout << endl
         << " Test case VII: j1 = " << j1 << ", j2 = " << j2
         << ", beta1 = " << beta1 << "/a0, beta2 = " << beta2
         << "/a0, R2-R1 = " << x0 << " a0";
    omegachi1 = pow(hydroconst::hbarc_eV_nm * beta1 / hydroconst::a0_nm, 2) /
                hydroconst::muc2_eV;
    omegachi2 = pow(hydroconst::hbarc_eV_nm * beta2 / hydroconst::a0_nm, 2) /
                hydroconst::muc2_eV;
    omega1 = (2.0 * j1 + 1.0) * omegachi1;
    omega2 = (2.0 * j2 + 1.0) * omegachi2;
    Areduced = 0.5 *
               pow(beta1 * hydroconst::hbarc_eV_nm / hydroconst::a0_nm, 2) /
               omegachi1 / hydroconst::muc2_eV;
    x0 *= hydroconst::a0_nm;
  } else {
    cout << endl
         << " Give the atomic mass of atom 1 (u) ....................: ";
    cin >> A1;
    cout << endl
         << " Give the atomic mass of atom 2 (u) ....................: ";
    cin >> A2;
    Areduced = A1 * A2 / (A1 + A2);
    cout << endl
         << " Give equilibrium distance of the lower oscillator (nm) : ";
    cin >> Re1;
    cout << endl
         << " Give hbar*omega of the lower oscillator (eV) ..........: ";
    cin >> omega1;
    cout << endl
         << " Give hbar*omegaChi of the lower oscillator (eV) .......: ";
    cin >> omegachi1;
    cout << endl
         << " Give equilibrium distance of the upper oscillator (nm) : ";
    cin >> Re2;
    cout << endl
         << " Give hbar*omega of the upper oscillator (eV) ..........: ";
    cin >> omega2;
    cout << endl
         << " Give hbar*omegaChi of the upper oscillator (eV) .......: ";
    cin >> omegachi2;

    x0 = Re2 - Re1;

    double r1 = 0.5 * omega1 / omegachi1;
    int vmax1 = int(r1 * (sqrt(2.0) + 1) - 0.49);
    double r2 = 0.5 * omega2 / omegachi2;
    int vmax2 = int(r2 * (sqrt(2.0) + 1) - 0.49);

    cout << endl << " Give max lower vibrational quantum number (<= ";
    cout.width(3);
    cout << vmax1 << ") ....: ";
    cin >> v1max;
    cout << endl << " Give max upper vibrational quantum number (<= ";
    cout.width(3);
    cout << vmax2 << ") ....: ";
    cin >> v2max;
    cout << endl;
  }

  int choice;
  cout << endl << " Options for additional output: ";
  cout << endl << "  0: no additional output";
  cout << endl << "  1: additional tabulation of transition energies";
  cout << endl << "  2: test of Schmidt recursion";
  cout << endl << "  3: harmonic FC factors according to Frank et al.";
  cout << endl << "  4: a=1 test of Morse integral";
  cout << endl << " Make a choice (0/1/2/3/4) .............................: ";
  cin >> choice;
  bool Eout_flag = (choice == 1);
  bool test_Schmidt_flag = (choice == 2);
  bool test_Frank_flag = (choice == 3);
  bool a1_test_flag = (choice == 4);

  cout << endl;

  double mreduced = hydroconst::muc2_eV * Areduced; // now in eV
  double mu = mreduced * pow(hydroconst::hbarc_eV_nm, -2);

  int v1dim = v1max + 1, v2dim = v2max + 1;
  matrix<double> Eharmo(v1dim, v2dim, 0.0);
  matrix<double> FCharmo(v1dim, v2dim, 0.0);
  matrix<double> Emorse(v1dim, v2dim, 0.0);
  matrix<double> FCmorse(v1dim, v2dim, 0.0);

  //--------------------------------------------------------------------------------
  // harmonic potentials

  // overlap factors for two displaced (by distance x0) harmonic oscillators
  // according to P. P. Schmidt, Mol. Phys. 108 (2010) 1513
  double a1 = mu * omega1;
  double a2 = mu * omega2;

  double ama = a1 - a2;
  double apa = a1 + a2;
  double fac1 = -a1 * sqrt(0.5 * a2) * x0 / apa;
  double fac2 = -ama / apa;
  double fac3 = sqrt(0.5 * a1) * a2 * x0 / apa;
  double fac4 = 2.0 * sqrt(a1 * a2) / apa;
  double fac5 = -fac3;
  double fac6 = -fac2;
  double fac7 = -fac1;

  matrix<double> I(v1dim, v2dim, 0.0);
  I(0, 0) = sqrt(2.0 * sqrt(a1 * a2) / apa) *
            exp(-0.5 * a1 * a2 * x0 * x0 / apa); // I_0,0, Eq. 19 of Schmidt
  if (v1max > 0) {
    I(1, 0) =
        sqrt(2.0 * a1) * a2 / apa * x0 * I(0, 0); // I_1,0, Eq. 20 of Schmidt
    for (int v1 = 0; v1 < v1max - 1; v1++) {
      // I_v1+2,0, Eq. 17 of Schmidt
      I(v1 + 2, 0) =
          sqrt((v1 + 1.0) / (v1 + 2.0)) * ama / apa * I(v1, 0) +
          sqrt(2.0 / (v1 + 2.0)) * sqrt(a1) * a2 / apa * x0 * I(v1 + 1, 0);
    }
  }
  if (v2max > 0) {
    I(0, 1) =
        -a1 * sqrt(2.0 * a2) / apa * x0 * I(0, 0); // I_0,1, Eq. 21 of Schmidt
    for (int v2 = 0; v2 < v2max - 1; v2++) {
      // I_0,v2+2, Eq. 18 of Schmidt
      I(0, v2 + 2) =
          -sqrt((v2 + 1.0) / (v2 + 2.0)) * ama / apa * I(0, v2) -
          sqrt(2.0 / (v2 + 2.0)) * a1 * sqrt(a2) / apa * x0 * I(0, v2 + 1);
    }
  }
  int nmax = (v1max > v2max) ? v1max : v2max;
  for (int n = 0; n < nmax; n++) {
    // recursion according to Eq. 12 of Schmidt
    for (int v1 = n; (n < v2max) && (v1 < v1max); v1++) {
      int v2 = n;
      if (test_Schmidt_flag) {
        cout << n << "    (" << v1 + 1 << "," << v2 + 1 << "):   (" << v1 + 1
             << "," << v2 << ")   (" << v1 + 1 << "," << v2 - 1 << ")   (" << v1
             << "," << v2 + 1 << ")   (" << v1 << "," << v2 << ")   (" << v1
             << "," << v2 - 1 << ")   (" << v1 - 1 << "," << v2 + 1 << ")   ("
             << v1 - 1 << "," << v2 << ")   (" << v1 - 1 << "," << v2 - 1
             << ")    ";
      }
      double s1 = 1.0 / sqrt(v1 + 1.0);
      double s2 = 1.0 / sqrt(v2 + 1.0);
      double nextI = fac1 * s2 * I(v1 + 1, v2);
      nextI += (v2 == 0) ? 0.0 : fac2 * sqrt(1.0 * v2) * s2 * I(v1 + 1, v2 - 1);
      nextI += fac3 * s1 * I(v1, v2 + 1);
      nextI += fac4 * s1 * s2 * I(v1, v2);
      nextI += v2 == 0 ? 0.0 : fac5 * sqrt(1.0 * v2) * s1 * s2 * I(v1, v2 - 1);
      nextI += v1 == 0 ? 0.0 : fac6 * sqrt(1.0 * v1) * s1 * I(v1 - 1, v2 + 1);
      nextI += v1 == 0 ? 0.0 : fac7 * sqrt(1.0 * v1) * s1 * s2 * I(v1 - 1, v2);
      nextI += ((v1 == 0) || (v2 == 0))
                   ? 0.0
                   : sqrt(1.0 * v1 * v2) * s1 * s2 * I(v1 - 1, v2 - 1);
      I(v1 + 1, v2 + 1) = nextI;
      if (test_Schmidt_flag)
        cout << nextI << endl;
    } // end for (v1...)
    for (int v2 = n + 1; (n < v1max) && (v2 < v2max); v2++) {
      int v1 = n;
      if (test_Schmidt_flag) {
        cout << n << "    (" << v1 + 1 << "," << v2 + 1 << "):   (" << v1 + 1
             << "," << v2 << ")   (" << v1 + 1 << "," << v2 - 1 << ")   (" << v1
             << "," << v2 + 1 << ")   (" << v1 << "," << v2 << ")   (" << v1
             << "," << v2 - 1 << ")   (" << v1 - 1 << "," << v2 + 1 << ")   ("
             << v1 - 1 << "," << v2 << ")   (" << v1 - 1 << "," << v2 - 1
             << ")    ";
      }
      double s1 = 1.0 / sqrt(v1 + 1.0);
      double s2 = 1.0 / sqrt(v2 + 1.0);
      double nextI = fac1 * s2 * I(v1 + 1, v2);
      nextI += (v2 == 0) ? 0.0 : fac2 * sqrt(1.0 * v2) * s2 * I(v1 + 1, v2 - 1);
      nextI += fac3 * s1 * I(v1, v2 + 1);
      nextI += fac4 * s1 * s2 * I(v1, v2);
      nextI += v2 == 0 ? 0.0 : fac5 * sqrt(1.0 * v2) * s1 * s2 * I(v1, v2 - 1);
      nextI += v1 == 0 ? 0.0 : fac6 * sqrt(1.0 * v1) * s1 * I(v1 - 1, v2 + 1);
      nextI += v1 == 0 ? 0.0 : fac7 * sqrt(1.0 * v1) * s1 * s2 * I(v1 - 1, v2);
      nextI += ((v1 == 0) || (v2 == 0))
                   ? 0.0
                   : sqrt(1.0 * v1 * v2) * s1 * s2 * I(v1 - 1, v2 - 1);
      I(v1 + 1, v2 + 1) = nextI;
      if (test_Schmidt_flag)
        cout << nextI << endl;
    } // end for(v2 ...)
  } // end for (n...)
  double E00 = 0.5 * (omega2 - omega1);
  cout << endl
       << " Franck-Condon factors for harmonic potentials according to "
          "Schmidt, Mol. Phys 108 (2010) 1513";
  cout << endl
       << " -------------------------------------------------------------------"
          "--------------------------";
  cout << endl << " v1\\v2 ";
  for (int v2 = 0; v2 <= v2max; v2++) {
    cout << "  ";
    cout.width(9);
    cout << v2;
  }
  cout << endl;
  double Emin = 0.0, Emax = 0.0;
  for (int v1 = 0; v1 <= v1max; v1++) {
    cout << " ";
    cout.width(5);
    cout << v1 << ":  FCf ";
    // if (fabs(omega1-omega2) < 1E-9) cout <<
    // pow(b1*b2,v2)*exp(-b1*b2-lnf(v2)); // Poisson distribution for v1=0 if
    // omega1=omega
    for (int v2 = 0; v2 <= v2max; v2++) {
      double fcf = pow(I(v1, v2), 2);
      cout << "  ";
      cout.precision(6);
      cout.width(9);
      cout << right << fixed << fcf;
      FCharmo(v1, v2) = fcf;
      FCmorse(v1, v2) = fcf; // default may be changed below
    }
    if (Eout_flag)
      cout << endl << "       DeltaE";
    for (int v2 = 0; v2 <= v2max; v2++) {
      if (Eout_flag) {
        cout << "  ";
        cout.precision(5);
        cout.width(9);
      }
      double E21 = (v2 + 0.5) * omega2 - (v1 + 0.5) * omega1;
      if (Eout_flag)
        cout << right << fixed << E21 - E00;
      Eharmo(v1, v2) = E21;
      Emorse(v1, v2) = E21; // default, may be changed below
      if (E21 > Emax)
        Emax = E21;
      if (E21 < Emin)
        Emin = E21;
    }
    cout << endl;
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // second method for comparison
  // overlap factors for two displaced (by distance x0) harmonic oscillators
  // according to Frank et al., J. Math. Chem. 25 (1999) 383 (numerically
  // unstable for high n)
  if (test_Frank_flag) {
    double so = omega1 + omega2;
    double po = omega1 * omega2;
    double sqrtpo = sqrt(po);

    // definitions according to Eqs. 10 & 12 of Frank et al.
    double u = 0.5 * so / sqrtpo; //  (omega2+omega1)/(2*sqrt(omega1*omega2))
    double v = 0.5 * (omega2 - omega1) /
               sqrtpo; //  (omega2-omega1)/(2*sqrt(omega1*omega2))
    double b1 = sqrt(0.5 * mu * omega1) * x0;
    double b2 = sqrt(0.5 * mu * omega2) * x0;

    // definitions according to Eq. 23 of Frank et al.
    double g = v / u;
    double f = (v * b2 + b1) / (u * u);
    double A00 =
        exp(-0.5 * mu * po * x0 * x0 / so) / sqrt(u); // Eq. 4 of Frank et al.

    cout << endl
         << " Franck-Condon factors for harmonic potentials according to Frank "
            "et al., J. Math. Chem. 25 (1999) 383 "
            "(numerically unstable for high v) ";
    cout << endl
         << " -----------------------------------------------------------------"
            "----------------------------------------"
            "-------------------------------";
    cout << endl << " v1\\v2 ";
    for (int v2 = 0; v2 <= v2max; v2++) {
      cout << "  ";
      cout.width(9);
      cout << v2;
    }
    cout << endl;
    for (int v1 = 0; v1 <= v1max; v1++) {
      cout << " ";
      cout.width(5);
      cout << v1 << ":  FCf ";
      // if (fabs(omega1-omega2) < 1E-9) cout <<
      // pow(b1*b2,v2)*exp(-b1*b2-lnf(v2)); // Poisson distribution for v1=0 if
      // omega1=omega
      for (int v2 = 0; v2 <= v2max; v2++) {
        double fcf = 0.0;
        double fac0 = sqrt(exp(lgamma(v2 + 1.0) - lgamma(v1 + 1.0)));
        for (int k = 0; k <= v2; k++) {
          int v2k = v2 - k;
          double powb2v2k = pow(b2, v2k);
          for (int s = 0; s <= k / 2; s++) {
            int sign = pow(-1, v2k + s);
            for (int t = 0; t <= (k - 2 * s); t++) {
              int st = s + t;
              int kst = k - st;
              int k2st = kst - s;
              double fac1k = powb2v2k * pow(u, kst) * pow(v, st) * pow(0.5, s);
              double fac2k = exp(-lgamma(v2k + 1.0) - lgamma(s + 1.0) -
                                lgamma(k2st + 1.0) - lgamma(t + 1.0));
              double sum = 0.0;
              int jmax = (v1 + t - k2st) / 2;
              for (int j = 0; j <= jmax; j++) {
                // Eq. 25 of A. Frank et al., J. Math. Chem. 25 (1999) 383
                int kk = v1 + t - k2st - 2 * j;
                sum += pow(f, kk) * pow(-0.5 * g, j) *
                       exp(lgamma(v1 + t + 1.0) - lgamma(kk + 1.0) -
                           lgamma(j + 1.0));
              }
               fcf += sign * fac0 * fac1k * fac2k * sum *
                     A00; // harmonicFC_A(v1+t,k2st,f,g,A00);
            } // end for(t ...)
          } // end for(s ...)
        } // end for(k ...)
        cout << "  ";
        cout.precision(6);
        cout.width(9);
        cout << right << fixed << fcf * fcf;
      }
      cout << endl;
    }
  } // end of harmonic calculation according to Frank et al.
  //--------------------------------------------------------------------------------
  if ((omegachi1 > 0.0) && (omegachi2 > 0.0)) {
    FCsetup_integration();
    cout << endl
         << " Franck-Condon factors for Morse potentials according to Lopez et "
            "al., Int. J. Quant. Chem, 88 (2002) 280";
#ifdef useBOOST
    cout << " (using BOOST extended precision arithmetic)";
#endif
    cout << endl
         << " -----------------------------------------------------------------"
            "---------------------------------------";
#ifdef useBOOST
    cout << "--------------------------------------------";
#endif
    if (a1_test_flag)
      cout << endl << " Test run: Output values are relative errors" << endl;
    cout << endl << " v1\\v2 ";
    for (int v2 = 0; v2 <= v2max; v2++) {
      cout << "  ";
      cout.width(9);
      cout << v2;
    }
    cout << endl;

    E00 = 0.5 * (omega2 - 0.5 * omegachi2) - 0.5 * (omega1 - 0.5 * omegachi1);
    for (int v1 = 0; v1 <= v1max; v1++) {
      cout << " ";
      cout.width(5);
      cout << v1 << ":  FCf ";
      for (int v2 = 0; v2 <= v2max; v2++) {
        cout << "  ";
        cout.precision(6);
        cout.width(9);
        double fcf, fcf1;
        if (a1_test_flag) {
          fcf = pow(MorseFC(v1, v2, omega1, omega2, omegachi1, omegachi2,
                            mreduced, x0, 1),
                    2);
          fcf1 = pow(MorseFC(v1, v2, omega1, omega2, omegachi1, omegachi2,
                             mreduced, x0, 2),
                     2);
          cout << right << scientific << fcf / fcf1 - 1.0;
        } else {
          fcf = pow(MorseFC(v1, v2, omega1, omega2, omegachi1, omegachi2,
                            mreduced, x0, 0),
                    2);
          cout << right << fixed << fcf;
          FCmorse(v1, v2) = fcf;
        }
      }
      if (Eout_flag)
        cout << endl << "       DeltaE";
      for (int v2 = 0; v2 <= v2max; v2++) {
        if (Eout_flag) {
          cout << "  ";
          cout.precision(5);
          cout.width(9);
        }
        double E21 = (v2 + 0.5) * (omega2 - (v2 + 0.5) * omegachi2) -
                     (v1 + 0.5) * (omega1 - (v1 + 0.5) * omegachi1);
        if (Eout_flag)
          cout << right << fixed << E21 - E00;
        Emorse(v1, v2) = E21;
        if (E21 > Emax)
          Emax = E21;
        if (E21 < Emin)
          Emin = E21;
      }
      cout << endl;
    }
    if (selected_test > 0) {
      cout << endl;
      cout << " omega1 = " << omega1 << " eV, omegachi1 = " << omegachi1
           << " eV" << endl;
      cout << " omega2 = " << omega2 << " eV, omegachi2 = " << omegachi2
           << " eV" << endl;
      cout << " x0 = " << x0 << " nm" << endl;
    }
  } // end  if ((omegachi1 > 0.0) && (omegachi2 > 0.0))

  cout << endl << " Error analysis? (y/n) .................................:  ";
  cin >> answer;
  bool error_analysis_flag = ((answer == 'y') || (answer == 'Y'));
  double domega2 = 0.0, domegachi2 = 0.0, dx0 = 0.0;
  matrix<double> FCmorseMin(v1dim, v2dim, 0.0);
  matrix<double> FCmorseMax(v1dim, v2dim, 0.0);
  if (error_analysis_flag) {
    cout << endl
         << " Give uncertainty of hbar*omega upper (eV) .............:  ";
    cin >> domega2;
    cout << endl
         << " Give uncertainty of hbar*omegachi upper (eV) ..........:  ";
    cin >> domegachi2;
    cout << endl
         << " Give uncertainty of Re upper (nm) .....................:  ";
    cin >> dx0;

    for (int v1 = 0; v1 <= v1max; v1++) {
      for (int v2 = 0; v2 <= v2max; v2++) {
        FCmorseMin(v1, v2) = FCmorse(v1, v2);
        FCmorseMax(v1, v2) = FCmorse(v1, v2);
      }
    }
    double Eomega2, Eomegachi2, Ex0;
    cout << endl
         << endl
         << " Performing error estimation of Franck-Condon factors";
    for (Eomega2 = omega2 - domega2; Eomega2 < omega2 + 1.1 * domega2;
         Eomega2 += domega2) {
      for (Eomegachi2 = omegachi2 - domegachi2;
           Eomegachi2 < omegachi2 + 1.1 * domegachi2;
           Eomegachi2 += domegachi2) {
        for (Ex0 = x0 - dx0; Ex0 < x0 + 1.1 * dx0; Ex0 += dx0) {
          cout << endl
               << " omgea: " << Eomega2 << "    omgeachi: " << Eomegachi2
               << "    Re: " << Re1 + Ex0;
          for (int v1 = 0; v1 <= v1max; v1++) {
            for (int v2 = 0; v2 <= v2max; v2++) {
              double fcf;
              fcf = pow(MorseFC(v1, v2, omega1, Eomega2, omegachi1, Eomegachi2,
                                mreduced, Ex0, 0),
                        2);
              if (fcf > FCmorseMax(v1, v2))
                FCmorseMax(v1, v2) = fcf;
              if (fcf < FCmorseMin(v1, v2))
                FCmorseMin(v1, v2) = fcf;
            }
          }
        }
      }
    }
    cout << endl;
  } // end if (error_analysis_flag)
  /////////////////////////////////////////////////////////////////////////////////////////

  cout << endl << " Output of vibrational spectra? (y/n) ..................: ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    int Jmax = 0;     // max rotational quantum number
    double Be1 = 0.0; // rotational constant of lower oscillator
    int Lambda1 = 0,
        Lambda2 = 0; // electronic angular momenta of lower and upper level
    int nucl_symmetry =
        0; // symmetry of nuclear wavefunction, 0: symmetric, 1: antisymmetric
    int elec_symmetry = 0; // symmetry of electronic wavefunction, 0: symmetric,
                           // 1: antisymmetric
    double Inucl = 0;      // nuclear spin

    fstream fout;
    string filename, outfilename;

    //  --- time stamp ---
    time_t rawtime;
    struct tm *timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    cout << endl
         << " Give filename for output (*.fcf) ......................: ";
    cin >> filename;
    outfilename = filename + ".fcf";
    fout.open(outfilename, fstream::out);
    fout << "##################################################################"
            "###########################"
         << endl;
    fout << "### Vibrational spectra from Franck-Condon factors of the "
            "harmonic and the Morse potentials "
         << endl;
#ifdef useBOOST
    fout << "### used BOOST multiprecision arithmetic for Morse Franck-Condon "
            "factors"
         << endl;
#endif
    fout << "###  hydrocal revision     : " << HYDROCAL_REVISION
         << endl; // defined in buildinfo.h
    fout << "###               filename : " << outfilename << endl;
    fout << "###      start date & time : "
         << asctime(timeinfo); // endl comes with asctime
    fout << "###     mass of atom 1 (u) : " << A1 << endl;
    fout << "###     mass of atom 2 (u) : " << A2 << endl;
    fout << "###                 max v1 : " << v1max << endl;
    fout << "###                 max v2 : " << v2max << endl;
    double kT = 0.0;
    cout << endl
         << " Give temperature in K .................................: ";
    cin >> kT;
    fout << "###                 T  (K) : " << kT << endl;
    kT *= hydroconst::kB_eV_K; // conversion K -> eV
    fout << "###                kT (eV) : " << kT << endl;
    cout << endl
         << " Give max rotational quantum number (0 ignores rotions) : ";
    cin >> Jmax;
    if (Jmax > 0) {
      fout << "###   max rot. quantum no. : " << Jmax << endl;
      cout << endl
           << " Give rotational constant Be of lower oscillator (cm⁻1) : ";
      cin >> Be1;
      fout << "###             Be1 (1/cm) : " << Be1 << endl;
      fout << "###             Be2 (1/cm) : " << Be1 * Re1 * Re1 / Re2 / Re2
           << endl;
      cout << endl
           << " Give angular momentum LAMBDA1 of lower level ..........: ";
      cin >> Lambda1;
      fout << "###  ang. momentum LAMBDA1 : " << Lambda1 << endl;
      cout << endl
           << " Give angular momentum LAMBDA2 of upper level ..........: ";
      cin >> Lambda2;
      fout << "###  ang. momentum LAMBDA2 : " << Lambda2 << endl;
      cout << endl
           << " Give nuclear spin I ...................................: ";
      cin >> Inucl;
      fout << "###           nuclear spin : " << Inucl << endl;
      nucl_symmetry =
          ((int)(2 * Inucl + 0.1) %
           2); // 0 for integer nuclear spin, 1 for half-integer nuclear spin
      if (nucl_symmetry == 0) {
        fout << "###             nuclei are : bosons" << endl;
      } else {
        fout << "###             nuclei are : fermionons" << endl;
      }

      cout << endl
           << " electronic wave function (anti)symmetric ? (a/s) ......: ";
      cin >> answer;
      if ((answer == 'a') || (answer == 'A')) {
        elec_symmetry = 1;
        fout << "###   elec. wavefunction : antisymmetric" << endl;
      } else {
        elec_symmetry = 0;
        fout << "###   elec. wavefunction : symmetric" << endl;
      }
    } // end if (Jmax>0)

    vector<double> wharmo(v1dim);
    vector<double> wmorse(v1dim);
    if (v1max > 0) {
      cout << endl
           << " Boltzmann population also of vibrational levels? (y/n) : ";
      cin >> answer;
      bool maxwell_flag = ((answer == 'Y') || (answer == 'y'));
      if (maxwell_flag) {
        double wsumh = 0.0, wsumm = 0.0;
        for (int v1 = 0; v1 <= v1max; v1++) {
          double Ev1 = (v1 + 0.5) * omega1;
          wharmo[v1] = exp(-Ev1 / kT);
          wsumh += wharmo[v1];
          Ev1 -= pow(v1 + 0.5, 2) * omegachi1;
          wmorse[v1] = exp(-Ev1 / kT);
          wsumm += wmorse[v1];
        } // end for(v1 ...)
        for (int v1 = 0; v1 <= v1max; v1++) {
          wharmo[v1] /= wsumh;
          wmorse[v1] /= wsumm;
        }
      } // end if (maxwell_flag)
      else {
        // manual input of initial level populations
        double wsum = 1.0;
        for (int v1 = 1; v1 <= v1max; v1++) {
          double wini;
          cout << endl
               << " Give population of vibrational level " << v1
               << " ..............: ";
          cin >> wini;
          wsum -= wini;
          wharmo[v1] = wini;
          wmorse[v1] = wini;
        } // end for (v1... )
        if (wsum < 0.0) {
          cout << endl
               << " ERROR: Sum of initial populations > 1 " << endl
               << endl;
          return 0;
        }
        wharmo[0] = wsum;
        wmorse[0] = wsum;
      } // end else (maxwell_flag)
    } // end if (v1max > 0)
    else {
      // only v1=0 is populated
      wharmo[0] = 1.0;
      wmorse[0] = 1.0;
    }
    for (int v1 = 0; v1 <= v1max; v1++) {
      fout << "###  ini pop harmo level " << v1 << " : " << wharmo[v1] << endl;
    }
    for (int v1 = 0; v1 <= v1max; v1++) {
      fout << "###  ini pop Morse level " << v1 << " : " << wmorse[v1] << endl;
    }

    // rotational energy shifts and relative line strengths
    double Erot1 = Be1 * 1E-7 * hydroconst::hc_eV_nm; // conversion cm-1 -> eV
    double Erot2 = Erot1 * Re1 * Re1 / Re2 / Re2;
    double wrot_sum = 0.0;
    vector<double> wrot(Jmax + 1);
    vector<double> SrotP(Jmax + 1);
    vector<double> SrotQ(Jmax + 1);
    vector<double> SrotR(Jmax + 1);
    vector<double> DErotP(Jmax + 1);
    vector<double> DErotQ(Jmax + 1);
    vector<double> DErotR(Jmax + 1);
    for (int J = 0; (Jmax > 0) && (J <= Jmax);
         J++) { // executed only for Jmax>0
      double JJ1 = J * (J + 1.0);
      wrot[J] = Inucl + (nucl_symmetry + elec_symmetry + J + 1) %
                            2; // nuclear spin statistics
      wrot[J] *= (2.0 * J + 1.0) *
                 exp(-JJ1 * Erot1 / kT); // factor (2J+1) is the statistical
                                         // weight of the rotational level
      // rotational energy shifts
      DErotP[J] = (J - 1.0) * J * Erot2 - JJ1 * Erot1;         // J'=J-1
      DErotQ[J] = JJ1 * (Erot2 - Erot1);                       //  J'=J
      DErotR[J] = (J + 1.0) * (J + 2.0) * Erot2 - JJ1 * Erot1; // J' = J + 1
      // Hoenl-London factors [Hansson & Watson, J. Mol. Spec. 233, 169 (2005)]
      if (Lambda2 == Lambda1) { // Delta Lambda = 0
        SrotP[J] = (J > 0) ? 1.0 * (J + Lambda1) * (J - Lambda1) / J : 0.0;
        SrotQ[J] = (J > 0) ? (2.0 * J + 1) * Lambda2 * Lambda2 / JJ1 : 0.0;
        SrotR[J] = 1.0 * (J + 1 + Lambda2) * (J + 1 - Lambda2) / (J + 1);
      } else if (Lambda2 > Lambda1) { // Delta Lambda = +1
        SrotP[J] = (J > 0) ? 0.25 * (J - 1 - Lambda1) * (J - Lambda1) / J : 0.0;
        SrotQ[J] = (J > 0) ? 0.25 * (J + Lambda2) * (J + 1 - Lambda2) *
                                 (2.0 * J + 1.0) / JJ1
                           : 0.0;
        SrotR[J] = 0.25 * (J + 1 + Lambda2) * (J + Lambda2) / (J + 1);
      } else { // Delta Lambda = -1
        SrotP[J] = (J > 0) ? 0.25 * (J - 1 + Lambda1) * (J + Lambda1) / J : 0.0;
        SrotQ[J] = (J > 0) ? 0.25 * (J - Lambda2) * (J + 1.0 + Lambda2) *
                                 (2.0 * J + 1.0) / JJ1
                           : 0.0;
        SrotR[J] = 0.25 * (J + 1 - Lambda2) * (J - Lambda2) / (J + 1);
      }
      SrotP[J] *= wrot[J];
      SrotQ[J] *= wrot[J];
      SrotR[J] *= wrot[J];
      wrot_sum += (SrotP[J] + SrotQ[J] + SrotR[J]);
    }
    if (wrot_sum > 0.0) {
      for (int J = 0; J <= Jmax; J++) {
        SrotP[J] /= wrot_sum;
        SrotQ[J] /= wrot_sum;
        SrotR[J] /= wrot_sum;
      }
    }

    double wLorentz, wGauss, DeltaTe;
    cout << endl
         << " Give lorentzian width in eV ...........................: ";
    cin >> wLorentz;
    cout << endl
         << " Give gaussian width in eV .............................: ";
    cin >> wGauss;
    cout << endl
         << " Give energy difference between electronic states in eV : ";
    cin >> DeltaTe;

    fout << "###               Re1 (nm) : " << Re1 << endl;
    fout << "###            omega1 (eV) : " << omega1 << endl;
    fout << "###         omegachi1 (eV) : " << omegachi1 << endl;
    fout << "###            omega2 (eV) : " << omega2 << endl;
    fout << "###         omegachi2 (eV) : " << omegachi2 << endl;
    fout << "###               Re2 (nm) : " << Re2 << endl;
    fout << "###           Re2-Re1 (nm) : " << x0 << endl;
    fout << "###           Te2-Te1 (eV) : " << DeltaTe << endl;
    fout << "###          wLorentz (eV) : " << wLorentz << endl;
    fout << "###             wGauss(eV) : " << wGauss << endl;
    if (error_analysis_flag) {
      fout << "###      delta omega2 (eV) : " << domega2 << endl;
      fout << "###  delta  omegachi2 (eV) : " << domegachi2 << endl;
      fout << "###          delta R2 (nm) : " << dx0 << endl;
    }
    // parameters of the Morse potentials
    double MorseD1 = 0.25 * omega1 * omega1 / omegachi1;
    fout << "###          Morse D1 (eV) : " << MorseD1 << endl;

    double MorseBeta1 = 0.1 * sqrt(2.0 * mreduced * omegachi1) /
                        hydroconst::hbarc_eV_nm; // in A-1
    fout << "###      Morse beta1 (A-1) : " << MorseBeta1 << endl;

    double MorseD2 = 0.25 * omega2 * omega2 / omegachi2;
    fout << "###          Morse D2 (eV) : " << MorseD2 << endl;
    if (error_analysis_flag) {
      double deltaD2 = sqrt(pow(2.0 * MorseD2 / omega2 * domega2, 2) +
                            pow(MorseD2 / omegachi2 * domegachi2, 2));
      fout << "###          delta D2 (eV) : " << deltaD2 << endl;
    }

    double MorseBeta2 = 0.1 * sqrt(2.0 * mreduced * omegachi2) /
                        hydroconst::hbarc_eV_nm; // in A-1
    fout << "###     Morse beta2 (A-1) : " << MorseBeta2 << endl;
    if (error_analysis_flag) {
      double deltaBeta2 = 0.5 * MorseBeta2 / omegachi2 * domegachi2;
      fout << "###     delta beta2 (A-1) : " << deltaBeta2 << endl;
    }
    fout << "#################################################################"
         << endl;
    fout << "###  energy (eV) [1]  harmonic spectrum [2]   Morse spectrum [3]  "
            " Nrot [4]   DErotP [5]   SrotP[6]   "
            "DErotQ[7]   SrotQ [8]   DErotR [9]   SrotR [10]  v1 [11]   v2 "
            "[12]   Eharmo (eV) [13]   FCharmo[14]   "
            "Emorse (eV) [15]   FCmorse [16]";
    if (error_analysis_flag)
      fout << "   FCmoseNegErr [17]   FCmorsePosErr [18]";
    fout << endl;
    fout << "###--------------------------------------------------------------"
         << endl;

    cout << endl
         << " Differences between vibrational level energies range from "
         << DeltaTe + Emin << " to " << DeltaTe + Emax << " eV." << endl;

    double xlo, xhi, xdelta;
    cout << endl
         << " Give energy range in eV (Emin Emax Edelta) ............: ";
    cin >> xlo >> xhi >> xdelta;

    int npts = int((xhi - xlo) / xdelta + 1.1);
    if (npts < (v1max + 1) * (v2max + 1))
      npts = (v1max + 1) * (v2max + 1);
    int vv1 = 0, vv2 = 0;
    for (int n = 0; n < npts; n++) {
      double x = xlo + n * xdelta;
      double yharmo = 0.0, ymorse = 0.0;
      for (int v2 = 0; v2 <= v2max; v2++) {
        for (int v1 = 0; v1 <= v1max; v1++) {
          if (Jmax == 0) {
            yharmo += wharmo[v1] * voigt(x, DeltaTe + Eharmo(v1, v2),
                                         FCharmo(v1, v2), wLorentz, wGauss);
            ymorse += wmorse[v1] * voigt(x, DeltaTe + Emorse(v1, v2),
                                         FCmorse(v1, v2), wLorentz, wGauss);
          } else { // rotational broadening if Jmax>0
            double yyh = 0.0, yym = 0.0;
            for (int J = 0; J <= Jmax; J++) {
              if (SrotP[J] > 0) {
                yyh += SrotP[J] * voigt(x, DeltaTe + Eharmo(v1, v2) + DErotP[J],
                                        FCharmo(v1, v2), wLorentz, wGauss);
                yym += SrotP[J] * voigt(x, DeltaTe + Emorse(v1, v2) + DErotP[J],
                                        FCmorse(v1, v2), wLorentz, wGauss);
              }
              if (SrotQ[J] > 0) {
                yyh += SrotQ[J] * voigt(x, DeltaTe + Eharmo(v1, v2) + DErotQ[J],
                                        FCharmo(v1, v2), wLorentz, wGauss);
                yym += SrotQ[J] * voigt(x, DeltaTe + Emorse(v1, v2) + DErotQ[J],
                                        FCmorse(v1, v2), wLorentz, wGauss);
              }
              if (SrotR[J] > 0) {
                yyh += SrotR[J] * voigt(x, DeltaTe + Eharmo(v1, v2) + DErotR[J],
                                        FCharmo(v1, v2), wLorentz, wGauss);
                yym += SrotR[J] * voigt(x, DeltaTe + Emorse(v1, v2) + DErotR[J],
                                        FCmorse(v1, v2), wLorentz, wGauss);
              }
            }
            yharmo += yyh * wharmo[v1];
            ymorse += yym * wmorse[v1];
          } // end else (Jmax==0)
        } // end for(v1...)
      } // end for(v2...)
      fout << x << ", " << yharmo << ", " << ymorse << ",  ";
      if ((Jmax > 0) && (n <= Jmax)) {
        fout << n << ",  " << DErotP[n] << ",  " << SrotP[n] << ",  "
             << DErotQ[n] << ",  " << SrotQ[n] << ",  " << DErotR[n] << ",  "
             << SrotR[n];
      }
      if ((vv1 <= v1max) && (vv2 <= v2max)) {
        if ((Jmax == 0) || (n > Jmax))
          fout
              << "  ,  ,  ,  ,  ,  ,  "; // skip columns containing no
                                         // information on rotational parameters
        fout << ", " << vv1 << ", " << vv2 << ",  "
             << DeltaTe + Eharmo(vv1, vv2) << ", " << FCharmo(vv1, vv2);
        fout << ", " << DeltaTe + Emorse(vv1, vv2) << ", " << FCmorse(vv1, vv2);
        if (error_analysis_flag)
          fout << ",  " << FCmorseMin(vv1, vv2) - FCmorse(vv1, vv2) << ",  "
               << FCmorseMax(vv1, vv2) - FCmorse(vv1, vv2);
        fout << endl;
      } else {
        fout << endl;
      }
      vv2++;
      if (vv2 >= v2dim) {
        vv1++;
        vv2 = 0;
      }
    }
    // there might still some rotational parameters left
    for (int n = npts; (Jmax > 0) && (n <= Jmax); n++) {
      fout << "  ,  ,  ,  " << n << ",  " << DErotP[n] << ",  " << SrotP[n]
           << ",  " << DErotQ[n] << ",  " << SrotQ[n] << ",  " << DErotR[n]
           << ",  " << SrotR[n] << ",  ,  ,  ,  ,  , ";
      if (error_analysis_flag)
        fout << ",  ,  ";
      fout << endl;
    }
    fout.close();
    cout << endl
         << " Vibrational spectra written to " << outfilename << "." << endl;
  } // end output of spectra

  /////////////////////////////////////////////////////////////////////////////////////////////////
  return 1;
}
