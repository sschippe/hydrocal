/**
 * @file peakfunctions.cxx
 * @brief Coding of peak functions
 * SPDX-License-Identifier: MIT
 *
 * @author Stefan Schippers
 * @date 2013-04-23
 * @version $Id $
 *
 */

#include "peakfunctions.h"
#include "Faddeeva_w.h"
#include "hydroconst.h"
#include "hydromath.h"
#include "readxsec.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

#undef useBOOST
#if __has_include(<boost/math/quadrature/gauss_kronrod.hpp>)
#if __has_include(<boost/math/quadrature/tanh_sinh.hpp>)
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#define useBOOST
#endif
#endif

using namespace std;
using hydroconst::pi;

void sum_of_peaks(void) {
  int read_mode = 0, peak_shape = PEAK_undefined, npeaks = 0,
      number_of_levels = 1, level_number = 1;
  string peak_filename;

  printf("\n Create a spectrum as a sum of peaks.\n");
  npeaks = open_peak_file(peak_filename, read_mode, number_of_levels);
  if (npeaks < 0) {
    printf("\n File %s not found.\n", peak_filename.c_str());
    return;
  }
  printf("\n %d peak data sets detected\n", npeaks);

  vector<double> energy(npeaks, 0.0);
  vector<double> strength(npeaks, 0.0);
  vector<double> qFano(npeaks, 0.0);
  vector<double> wLorentz(npeaks, 0.0);
  vector<double> wGauss(npeaks, 0.0);
  double emin, emax;

  peak_shape =
      read_peak_file(peak_filename, read_mode, level_number, energy, strength,
                     qFano, wLorentz, wGauss, npeaks, emin, emax);
  if (peak_shape == PEAK_undefined)
    return;

  if (number_of_levels > 1) {
    printf("\n Data for %d different initial levels found.", number_of_levels);
    printf("\n Give the number (in the range 1 - %d) of the level to be used: ",
           number_of_levels);
    scanf("%d", &level_number);
    if ((level_number < 1) || (level_number > number_of_levels)) {
      printf("\n\n ERROR: Level number of of range!\n\n");
      exit(0);
    }
  }

  double wConvolution = 0.0;
  const double eps = 1E-9;
  if (peak_shape == PEAK_delta) {
    printf("\n Convolute delta peaks with Gaussian?");
    printf("\n  (no convolution for zero width)");
    printf("\n Give Gaussian width ..........: ");
    scanf("%lf", &wConvolution);
    if (wConvolution > 0.0)
      peak_shape = PEAK_Gauss;
  } else if ((peak_shape == PEAK_Lorentz) ||
             (peak_shape == PEAK_Lorentz_Steih)) {
    printf("\n Convolute Lorentzian peaks with Gaussian or P04 energy "
           "distribution?");
    printf("\n  (no convolution for zero width)");
    printf("\n Give width (>0 Gauss, <0 P04 exit slit in mum) ..........: ");
    scanf("%lf", &wConvolution);
    if (wConvolution > eps) {
      peak_shape = PEAK_Voigt;
    } else if (wConvolution < -eps) {
      peak_shape = PEAK_P04_Lorentz;
    }
  } else if (peak_shape == PEAK_Fano) {
    printf("\n Convolute Fano peaks with Gaussian?");
    printf("\n  (no convolution for zero width)");
    printf("\n Give Gaussian width ..........: ");
    scanf("%lf", &wConvolution);
    if (wConvolution > eps)
      peak_shape = PEAK_FanoVoigt;
  } else {
    printf("\n Convolute peaks with Gaussian?");
    printf("\n  (no convolution for zero width)");
    printf("\n Give Gaussian width ..........: ");
    scanf("%lf", &wConvolution);
  }
  if (wConvolution > eps) {
    for (int i = 0; i < npeaks; i++)
      wGauss[i] = sqrt(wGauss[i] * wGauss[i] + wConvolution * wConvolution);
  } else if (wConvolution < -eps) {
    for (int i = 0; i < npeaks; i++)
      wGauss[i] = -wConvolution;
  }

  double xx, xdelta, xmin, xmax;
  double yy;
  ofstream fout;
  string filenameout;

  printf("\n Peak positions range from %g to %g\n", emin, emax);
  printf("\n Give emin, emax, and edelta .............................: ");
  cin >> xmin >> xmax >> xdelta;
  int npts = int(fabs(xmax - xmin) / xdelta + 1.1);

  printf("\n Give name of output file ................................: ");
  cin >> filenameout;
  fout.open(filenameout);

  if (peak_shape == PEAK_Gauss) {
    printf("\n Calculating sum of Gauss peaks ...");
    for (int n = 0; n < npts; n++) {
      xx = xmin + n * xdelta;
      yy = 0.0;
      for (int i = 0; i < npeaks; i++)
        yy += gauss(xx, energy[i], strength[i], wGauss[i]);
      fout << setw(20) << xx << " " << setw(20) << yy << "\n";
    }
  } else if ((peak_shape == PEAK_Lorentz) ||
             (peak_shape == PEAK_Lorentz_Steih)) {
    printf("\n Calculating sum of Lorentz peaks ...");
    for (int n = 0; n < npts; n++) {
      xx = xmin + n * xdelta;
      yy = 0.0;
      for (int i = 0; i < npeaks; i++)
        yy += lorentz(xx, energy[i], strength[i], wLorentz[i]);
      fout << setw(20) << xx << " " << setw(20) << yy << "\n";
    }
  } else if (peak_shape == PEAK_Voigt) {
    printf("\n Calculating sum of Voigt peaks ...");
    for (int n = 0; n < npts; n++) {
      xx = xmin + n * xdelta;
      yy = 0.0;
      for (int i = 0; i < npeaks; i++)
        yy += voigt(xx, energy[i], strength[i], wLorentz[i], wGauss[i]);
      fout << setw(20) << xx << " " << setw(20) << yy << "\n";
    }
  } else if (peak_shape == PEAK_Fano) {
    printf("\n Calculating sum of Fano peaks ...");
    for (int n = 0; n < npts; n++) {
      xx = xmin + n * xdelta;
      yy = 0.0;
      for (int i = 0; i < npeaks; i++)
        yy += fano(xx, energy[i], strength[i], qFano[i], wLorentz[i]);
      fout << setw(20) << xx << " " << setw(20) << (double)yy << "\n";
    }
  } else if (peak_shape == PEAK_FanoVoigt) {
    printf("\n Calculating sum of FanoVoigt peaks ...");
    for (int n = 0; n < npts; n++) {
      xx = xmin + n * xdelta;
      yy = 0.0;
      for (int i = 0; i < npeaks; i++)
        yy += fanovoigt(xx, energy[i], strength[i], qFano[i], wLorentz[i],
                        wGauss[i]);
      fout << setw(20) << xx << " " << setw(20) << yy << "\n";
    }
  } else if (peak_shape == PEAK_P04_Lorentz) {
#ifdef useBOOST
    printf("\n Calculating sum of P04lorentz peaks ...");
    for (int n = 0; n < npts; n++) {
      double maxerror = 0.0, pemin = 0.0, pemax = 0.0;
      xx = xmin + n * xdelta;
      yy = 0.0;
      for (int i = 0; i < npeaks; i++) {
        double error;
        if (wLorentz[i] > 1.0E-9) {
          yy += P04lorentz(xx, energy[i], strength[i], wLorentz[i], wGauss[i],
                           error, pemin,
                           pemax); // wGauss[i] contains the exitslit width
          if (error > maxerror)
            maxerror = error;
        } else {
          yy += P04ped(xx, energy[i],
                       wGauss[i]); // wGauss[i] contains the exitslit width
        }
      }
      fout << setw(20) << xx << " " << setw(20) << yy << " " << setw(20)
           << maxerror << "\n";
    }
#else
    printf("\n Need BOOST library for numerical integration for P04lorentz "
           "peaks. If installed use \"make hydrocal-boost\"");
#endif
  }

  fout.close();
}

////////////////////////////////////////////////////
/**
 * @brief Area normalized gaussian
 *
 *	@param E  energy axis
 *  @param E0 peak center energy
 *  @param wg gaussian linewidths (FWHM)
 *	@param a  peak area
 *
 *  @return profile height at the given energy
 */
double gauss(double E, double E0, double a, double wg) {
  const double sqrtln2 = hydroconst::sqrtln2;
  const double sqrtpi = hydroconst::sqrtpi;
  double w = 0.5 * wg / sqrtln2;
  double eps = (E - E0) / w;
  if (fabs(eps) > 10.0)
    return 0.0;
  double eps2 = eps * eps;

  return a / w * exp(-eps2) / sqrtpi;
}

////////////////////////////////////////////////////
/**
 * @brief Area normalized lorentzian
 *
 *  @param E  energy axis
 *  @param E0 peak center energy
 *  @param a  peak area
 *  @param wl lorentzian linewidths (FWHM)
 *
 *  @return profile height at the given energy
 */
double lorentz(double E, double E0, double a, double wl) {

  double eps = 2.0 * (E - E0) / wl;

  return 2.0 * a / (1.0 + eps * eps) / wl / pi;
}

////////////////////////////////////////////////////
/**
 * @brief Fano profile
 *
 * @param E   energy axis
 * @param E0  peak center energy
 * @param a   peak area
 * @param q   esymmetry parameter
 * @param wl  lorentzian linewidths (FWHM)
 *
 * @return height of the profile at the given energy (goes to zero for
 * E->+-infinity)
 */
double fano(double E, double E0, double a, double q, double wl) {
  double eps = 2.0 * (E - E0) / wl;
  double fano = 2.0 * a / wl / pi *
                ((q + eps) * (q + eps) / (1.0 + eps * eps) - 1.0) /
                (q * q - 1.0);
  return fano;
}

/////////////////////////////////////////////////////////////////////
/**
 * @brief Voigt profile (Convolution of a lorentzian with a gaussian)
 *
 * @param E   energy axis
 * @param E0  peak center energy
 * @param a   peak area
 * @param wl  lorentzian linewidths (FWMH)
 * @param wg  gaussian linewidths (FWHM)
 *
 * @return height of the profile at the given energy (goes to zero for
 * E->+-infinity)
 */
double voigt(double E, double E0, double a, double wl, double wg) {

  if (fabs(wg) < 1e-20) {
    return lorentz(E, E0, a, wl);
  } else {
    const double sqrtln2 = hydroconst::sqrtln2;
    const double sqrtpi = hydroconst::sqrtpi;
    double x, y;
    x = 2.0 * sqrtln2 * (E - E0) / wg;
    y = sqrtln2 * wl / wg;

    // complex_error_function(x,y,re,im);
    std::complex<double> z(x, y);
    std::complex<double> cerf = Faddeeva::w(z);
    return a * real(cerf) * 2.0 * sqrtln2 / sqrtpi / wg;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Convolution of a Fano profile with a gaussian
 *
 * For details see Schippers, J. Quant. Spectrosc. Radiat. Transfer 219, 33 (2018);
 * https://doi.org/10.1016/j.jqsrt.2018.08.003.
 *
 * @param Fano profile convoluted with a gaussian
 * @param E  Energy axis
 * @param E0 Peak center energy
 * @param a  Peak area
 * @param q  Asymmetry parameter
 * @param wl Lorentzian linewidths (FWMH)
 * @param wg Gaussian linewidths (FWHM)
 *
 * @return Height of the profile at the given energy (goes to zero for
 * E->+-infinity)
 */
double fanovoigt(double E, double E0, double a, double q, double wl,
                 double wg) {
  double q2 = q * q;
  if (fabs(q) < 1e-10)
    q = 1e-10;
  if (q2 >= 1.0e6) { // symmetric peak
    return voigt(E, E0, a, wl, wg);
  } else { // asymmetric peak
    if (fabs(wg) < 1e-20) {
      return fano(E, E0, a, q, wl);
    } else {
      const double sqrtln2 = hydroconst::sqrtln2;
      const double sqrtpi = hydroconst::sqrtpi;

      double x, y, re, im;

      x = 2.0 * sqrtln2 * (E0 - E) / wg; // the sign matters!
      y = sqrtln2 * wl / wg;
      // complex_error_function(x,y,re,im);
      std::complex<double> z(x, y);
      std::complex<double> cerf = Faddeeva::w(z);
      re = real(cerf);
      im = imag(cerf);

      double qfac1 = q2 - 1.0;
      double qfac2 = 2.0 * q;
      double fac = 2.0 * fabs(a / qfac1) / sqrtpi / wl;

      return fac * y * (re * qfac1 - im * qfac2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Convolution of a Fano profile with a gaussian (returning also partial
 * derivatives with respecte to parameters)
 *
 * 2018-07-10: Attention needs to be updated because of q� -> q�-1 substitution
 *
 * @param E  energy axis
 * @param E0 peak center energy
 * @param a  peak area (in the limit q -> infinity)
 * @param q  asymmetry parameter
 * @param wl lorentzian linewidths (FWMH)
 * @param wg gaussian linewidths (FWHM)
 * @param df_de0 on exit, partial derivative of the profile with respect to E0
 * @param df_da  on exit, partial derivative of the profile with respect to a
 * @param df_dq  on exit, partial derivative of the profile with respect to q
 * @param df_dwl on exit, partial derivative of the profile with respect to wl
 * @param df_dwg on exit, partial derivative of the profile with respect to wg
 * @return profile height at the given energy (goes to zero for E->+-infinity)
 */
double fanovoigtderiv(double E, double E0, double a, double q, double wl,
                      double wg, double &df_dE0, double &df_da, double &df_dq,
                      double &df_dwl, double &df_dwg) {
  const double sqrtln2 = hydroconst::sqrtln2;
  const double sqrtpi = hydroconst::sqrtpi;

  if (fabs(q) < 1e-10)
    q = 1e-10;

  double q2 = q * q;
  if (q2 >= 1.0e6) {        // symmetric peak
    if (fabs(wg) < 1e-20) { // Lorentz
      double x = 2.0 * (E - E0) / wl;
      double fac = 2.0 * a / wl / pi;
      double f = fac / (1.0 + x * x);

      double df_dx = -2.0 * x * f / (1.0 + x * x);

      df_dE0 = -df_dx * 2.0 / wl;
      df_dq = 0.0;
      df_dwg = 0.0;
      df_dwl = -(f + x * df_dx) / wl;
      df_da = f / a;

      return f;
    } else { // Voigt
      double fac, x, y, re, im;

      fac = a * 2.0 * sqrtln2 / sqrtpi / wg;
      x = 2.0 * sqrtln2 * (E - E0) / wg;
      y = sqrtln2 * wl / wg;

      //	complex_error_function(x,y,re,im);
      std::complex<double> z(x, y);
      std::complex<double> cerf = Faddeeva::w(z);
      re = real(cerf);
      im = imag(cerf);

      double dre_dx = -2 * (x * re - y * im);
      double dre_dy = 2 * (y * re + x * im - 1 / sqrtpi);

      double f = fac * re;

      double df_dx = fac * dre_dx;
      double df_dy = fac * dre_dy;

      df_dE0 = -df_dx * 2.0 * sqrtln2 / wg;
      df_dq = 0.0;
      df_dwg = -(f + x * df_dx + y * df_dy) / wg;
      df_dwl = y * df_dy / wl;
      df_da = f / a;

      return f;
    }
  } else {                  // asymmetric peak
    if (fabs(wg) < 1e-20) { // Fano
      double x = 2.0 * (E - E0) / wl;
      double xq = 1.0 + x / q;
      double x2 = x * x + 1.0;
      double fac = 2.0 * a / wl / pi;
      double f = fac * (xq * xq / x2 - 1);

      double df_dx = 2 * xq / (x2 * x2) * (1.0 / q - x);

      df_dE0 = -2.0 * df_dx / wl;
      df_dq = -2.0 * fac * xq / (q2 * x2);
      df_dwg = 0.0;
      df_dwl = -(f + x * df_dx) / wl;
      df_da = f / a;

      return f;
    } else { // FanoVoigt
      double fac, qfac1, qfac2, x, y, re, im;

      fac = a * 2.0 * sqrtln2 / sqrtpi / wg;
      qfac1 = 1.0 - 1.0 / q2;
      qfac2 = 2.0 / q;
      x = 2.0 * sqrtln2 * (E0 - E) / wg;
      y = sqrtln2 * wl / wg;

      // complex_error_function(x,y,re,im);
      std::complex<double> z(x, y);
      std::complex<double> cerf = Faddeeva::w(z);
      re = real(cerf);
      im = imag(cerf);

      double dre_dx = -2 * (x * re - y * im);
      double dre_dy = 2 * (y * re + x * im - 1 / sqrtpi);
      double dim_dx = -dre_dy;
      double dim_dy = dre_dx;

      double f = fac * (re * qfac1 - im * qfac2);

      double df_dre = fac * qfac1;
      double df_dim = -fac * qfac2;
      double df_dx = df_dre * dre_dx + df_dim * dim_dx;
      double df_dy = df_dre * dre_dy + df_dim * dim_dy;

      df_dE0 = df_dx * 2.0 * sqrtln2 / wg;
      df_dq = 2.0 * fac * (re / q - im) / q2;
      df_dwg = -(f + x * df_dx + y * df_dy) / wg;
      df_dwl = y * df_dy / wl;
      df_da = f / a;

      return f;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Convolution of a lorentzian with a normalized trapezoidal
 *
 * @param E  energy axis
 * @param E0 peak center energy
 * @param a peak area
 * @param wl lorentzian width (FWHM)
 * @param wb trapezoidal base width
 * @param wt trapezoidal top width
 *
 * @return profile height at the given energy
 */
double trapezlorentzian(double E, double E0, double a, double wl, double wb,
                        double wt) {

  if (wb < wt)
    wb = wt;

  double peak = 2.0 * a / pi / (wb + wt) *
                (atan(2.0 * (E + 0.5 * wt - E0) / wl) -
                 atan(2.0 * (E - 0.5 * wt - E0) / wl));

  if (wt < wb) {
    double wfac = 1.0 / (wb * wb - wt * wt);
    peak += -4.0 * a / pi * wfac * (E - E0 - 0.5 * wb) *
            (atan(2.0 * (E - 0.5 * wt - E0) / wl) -
             atan(2.0 * (E - 0.5 * wb - E0) / wl));
    peak += 4.0 * a / pi * wfac * (E - E0 + 0.5 * wb) *
            (atan(2.0 * (E + 0.5 * wb - E0) / wl) -
             atan(2.0 * (E + 0.5 * wt - E0) / wl));
    peak += wfac * a * wl / pi *
            log((pow(E - 0.5 * wt - E0, 2) + pow(0.5 * wl, 2)) /
                (pow(E - 0.5 * wb - E0, 2) + pow(0.5 * wl, 2)));
    peak += -wfac * a * wl / pi *
            log((pow(E + 0.5 * wb - E0, 2) + pow(0.5 * wl, 2)) /
                (pow(E + 0.5 * wt - E0, 2) + pow(0.5 * wl, 2)));
  }
  return peak;
}
////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Rydberg series of FanoGauss peaks
 *
 * Quantum defect formula for photoionization cross section,
 * [Tully et al., Astron. Astrophys. 211 (1989) 485].
 * The experimental photon energy spread is accounted for by convolution with a
 * gaussian.
 *
 * @param x independent variable
 * @param dx high-n cut-off parameter (should be of the order of <i>wg</i>)
 * @param wG gaussian width
 * @param zeff effective nuclear charge
 * @param slim energy of series limit
 * @param m0,m1   quantum defect parameters of series
 * @param q0,q1   Fano asymmetry parameters of series
 * @param w0,w1   natural line width parameters of series
 *
 * @return resonance cross section at energy <i>x</i>
 */
double RydbergSeries(double x, double dx, double wG, double zeff, double slim,
                     double m0, double m1, double q0, double q1, double w0,
                     double w1) {
  const double Ryd = hydroconst::Ryd_eV;
  const double sqrtln2 = hydroconst::sqrtln2;
  ;
  const double sqrtpi = hydroconst::sqrtpi;

  if (x <= 0.0)
    return 0.0;

  double z2Ryd = zeff * zeff * Ryd;
  double eps = (x - slim) / z2Ryd;

  // if the Rydberg resonances are too closely spaced the numerical
  // evaluation becomes error prone. Therefore high Rydberg resonances are
  // cut off if two neighboring resonance positions are smaller than dx
  double epscut = dx > 0 ? -exp(2.0 * log(dx / (2 * z2Ryd)) / 3.0) : 0;
  if (eps > epscut)
    return 1.0;

  double m = m0 + eps * m1; // quantum defect
  double q = q0 + eps * q1; // asymmetry parameter
  double w = w0 + eps * w1; // width parameter
  if (w < 1e-20) {
    w = 1e-20;
  }
  double t = tanh(pi * w);
  double xx = tan(pi * (1.0 / sqrt(-eps) + m));

  double fano;

  if (wG == 0.0) {
    fano = (xx < 1e10 * t) ? pow(xx + q * t, 2) / (t * t + xx * xx)
                           : pow(1.0 + q * t / xx, 2);
  } else { // convolution of Fano profiles with a Gaussian
    double Rez, Imz, Rew, Imw;

    double swG = sqrtln2 / wG;
    Rez = -xx * swG;
    Imz = t * swG;
    // complex_error_function(Rez,Imz,Rew,Imw);
    std::complex<double> z(Rez, Imz);
    std::complex<double> cerf = Faddeeva::w(z);
    Rew = real(cerf);
    Imw = imag(cerf);

    fano = Imz * ((q * q - 1.0) * Rew - 2.0 * q * Imw) * sqrtpi + 1.0;
  }
  return fano;
}

//////////////////////////////////////////////////////
/**
 * @brief P04 photon-energy distribution
 *
 * @param E photon energy in eV
 * @param Enom nominal photon energy in eV where the distribution is "centered"
 * @param exitslit width of the monochromator exit slit in �m
 *
 * @return if (E>=0) line profile
 * @return if (E=-1) parameter x0 in eV
 * @return if (E=-4) parameter x3 in eV
 */
double P04ped(double E, double Enom, double exitslit) {
  int iselect = int(E - 0.1);
  double p0, p1, p2, p3;

  p0 = 0.99998 + exitslit * (-2.66898E-7 +
                             exitslit * (1.01814E-11 - exitslit * 1.43344E-14));
  p1 = -7.24403E-9 +
       exitslit *
           (-5.73555E-10 + exitslit * (-1.17566E-12 + exitslit * 1.69091E-15));
  p2 = 4.90297E-12 +
       exitslit *
           (1.64895E-13 + exitslit * (4.02126E-16 - exitslit * 6.94231E-19));
  p3 = -1.47683E-15 +
       exitslit *
           (-3.30761E-17 + exitslit * (-4.05937E-20 + exitslit * 1.12127E-22));
  double x0 = Enom * (p0 + Enom * (p1 + Enom * (p2 + Enom * p3)));
  if (iselect == -1)
    return x0;

  p0 = 1.00002 + exitslit * (1.69063E-7 +
                             exitslit * (2.51377E-10 - exitslit * 2.57658E-13));
  p1 = -1.8254E-8 +
       exitslit *
           (7.84451E-10 + exitslit * (6.90654E-13 - exitslit * 1.16341E-15));
  p2 = 2.15165E-11 +
       exitslit *
           (-3.23559E-13 + exitslit * (-1.34184E-16 + exitslit * 3.91997E-19));
  p3 = -5.71195E-15 +
       exitslit *
           (6.22553E-17 + exitslit * (2.10840E-20 - exitslit * 7.59955E-23));
  double x3 = Enom * (p0 + Enom * (p1 + Enom * (p2 + Enom * p3)));
  if (iselect == -4)
    return x3;

  p0 = 1.00002 + exitslit * (-2.63299E-7 +
                             exitslit * (6.47586E-11 - exitslit * 3.87743E-14));
  p1 = -2.99741E-8 +
       exitslit *
           (-5.06536E-10 + exitslit * (-1.63302E-12 + exitslit * 2.00357E-15));
  p2 = 2.60117E-11 +
       exitslit *
           (5.90389E-14 + exitslit * (9.9821E-16 - exitslit * 1.12591E-18));
  p3 = -7.36701E-15 +
       exitslit *
           (6.44697E-18 + exitslit * (-2.47334E-19 + exitslit * 2.65161E-22));
  double x1 = Enom * (p0 + Enom * (p1 + Enom * (p2 + Enom * p3)));

  p0 = 0.99999 + exitslit * (1.50239E-7 +
                             exitslit * (2.76282E-10 - exitslit * 2.96455E-13));
  p1 = -5.89483E-9 +
       exitslit *
           (7.89369E-10 + exitslit * (7.40758E-13 - exitslit * 1.13647E-15));
  p2 = 1.01221E-11 +
       exitslit *
           (-3.14101E-13 + exitslit * (-2.14376E-16 + exitslit * 3.94898E-19));
  p3 = -3.42145E-15 +
       exitslit *
           (6.42251E-17 + exitslit * (3.19619E-20 - exitslit * 7.0249E-23));
  double x2 = Enom * (p0 + Enom * (p1 + Enom * (p2 + Enom * p3)));

  p0 = 0.84812 + exitslit * (1.25075E-4 - exitslit * 3.94039E-8);
  p1 = 2.62001E-4 + exitslit * (-8.21621E-7 + exitslit * 6.48614E-10);
  p2 = -1.81687E-7 + exitslit * (8.98529E-10 - exitslit * 7.65058E-13);
  p3 = 4.367E-11 + exitslit * (-2.57965E-13 + exitslit * 2.23583E-16);
  double y1 = p0 + Enom * (p1 + Enom * (p2 + Enom * p3));

  p0 = 0.75495 + exitslit * (3.78988E-4 - exitslit * 8.17095E-8);
  p1 = 4.40277E-4 + exitslit * (-1.77432E-6 + exitslit * 9.64659E-10);
  p2 = -2.59264E-7 + exitslit * (1.54723E-9 - exitslit * 1.07759E-12);
  p3 = 5.02021E-11 + exitslit * (-4.08952E-13 + exitslit * 3.06506E-16);
  double y2 = p0 + Enom * (p1 + Enom * (p2 + Enom * p3));

  double normalization = 2.0 / ((x2 - x0) * y1 + (x3 - x1) * y2);

  double b0 = y1 / (x1 - x0);
  double b1 = (y2 - y1) / (x2 - x1);
  double b2 = -y2 / (x3 - x2);
  double a0 = -b0 * x0;
  double a1 = y1 - b1 * x1;
  double a2 = y2 - b2 * x2;

  if ((E >= x0) && (E < x1)) {
    return normalization * (a0 + b0 * E);
  } else if ((E >= x1) && (E < x2)) {
    return normalization * (a1 + b1 * E);
  } else if ((E >= x2) && (E < x3)) {
    return normalization * (a2 + b2 * E);
  } else {
    return 0.0;
  }
}

double Ex0(double E, double exitslit) {
  const double eps = 1e-6;
  double Emin = 240, Emax = 2000, Enew;
  if (E < Emin)
    return Emin + 0.007 * exitslit;
  if (E > Emax)
    return Emax + 0.007 * exitslit;
  while (fabs(Emax - Emin) > eps) {
    Enew = 0.5 * (Emin + Emax);
    double x0 = P04ped(-1, Enew, exitslit);
    if (x0 < E) {
      Emin = Enew;
    } else {
      Emax = Enew;
    }
  }
  return Enew;
}

double Ex3(double E, double exitslit) {
  const double eps = 1e-6;
  double Emin = 240, Emax = 2000, Enew;
  if (E < Emin)
    return Emin - 0.007 * exitslit;
  if (E > Emax)
    return Emax - 0.007 * exitslit;
  while (fabs(Emax - Emin) > eps) {
    Enew = 0.5 * (Emin + Emax);
    double x3 = P04ped(-4, Enew, exitslit);
    if (x3 < E) {
      Emin = Enew;
    } else {
      Emax = Enew;
    }
  }
  return Enew;
}

////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Convolution of a Lorentzian with the P04 photon-energy distribution
 *
 * @param E photon energy in eV
 * @param Eres resonance energy in eV
 * @param S resonance stregth in Mb eV
 * @param wL Lorentzian width
 * @param exitslit width of the P04 monochormator exit slit in mum
 * @param error on exit: error estimate from numerical integration
 *
 * @return profile height at the given energy
 */
double P04lorentz(double E, double Eres, double S, double wL, double exitslit,
                  double &error, double &Eintmin, double &Eintmax) {
#ifdef useBOOST
  const double Pi = hydroconst::pi;

  // note that the integrand is nonzero only for x0(Eint) <= E <= x3(Eint)
  Eintmin = Ex3(E, exitslit); // x3(Eint) determines the lowest value for Eint
  Eintmax = Ex0(E, exitslit); // x0(Eint) determines the highest value for Eint
  // the integrand is defined as an anonymous function (also called lambda
  // function) see
  // https://en.wikipedia.org/wiki/Anonymous_function#C++_(since_C++11)
  auto integrand = [&](double Eint) -> double {
    double x = 2.0 * (Eint - Eres) / wL;
    return P04ped(E, Eint, exitslit) / (x * x + 1.0);
  };
  double result = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(
      integrand, Eintmin, Eintmax, 7, 1e-12, &error);
  // boost::math::quadrature::tanh_sinh<double> integrator;
  // double result =integrator.integrate(integrand,Eintmin,Eintmax);
  return result * S / Pi * 2.0 / wL;
#else
  return 0.0;
#endif
}
