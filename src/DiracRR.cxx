/**
 * @file DiracRR.cxx
 *
 * @brief Relativistic (Dirac) hydrogenic radiative-recombination crosscsections
 *
 * @par CREATION
 * @author Stefan Schippers
 * @date 2026
 *
 * @par VERSION
 * @verbatim
 * $Id: DiracRR.cxx 2119 2026-09-22 14:52:53Z iamp $
 // SPDX-License-Identifier: MIT
 * @endverbatim
 */

#include "IchiharaEichlerRRdata.h"
#include "hydroconst.h"
#include "hydromath.h"
#include "hypergeometric.h"
#include "clebsch.h"
#include "DiracRate.h"
#include "DiracAngle.h"
#include "osci.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#if __has_include(<boost/numeric/odeint.hpp>)
#include <boost/numeric/odeint.hpp>
#define HYDROCAL_HAVE_BOOST_ODEINT 1
#endif


using namespace std;

namespace {

// Fixed state ordering used by the data tables, keyed by (n, kappa).
struct NK {
  int n, kappa;
};
  
const NK kStates[] = {
    {1, -1},  // 1s1/2
    {2, -1},  // 2s1/2
    {2, 1},   // 2p1/2
    {2, -2},  // 2p3/2
    {3, -1},  // 3s1/2
    {3, 1},   // 3p1/2
    {3, -2},  // 3p3/2
    {3, 2},   // 3d3/2
    {3, -3},  // 3d5/2
};

// Per-state natural cubic spline in log10(Te) vs log10(sigma*Te) space,
// built once for every tabulated (Z, state) on first use.
struct Interp {
  int Z; // nuclear charge
  int state; // one of the states defined just above 
  vector<double> lx, ly, y2;
};

//////////////////////////////////////////////////////////////////////////
/**
 * @brief builds the per-(Z, state) cubic-spline interpolators
 *
 * The interpolators are constructed once, on first use, for every
 * tabulated Z and state.  Each interpolates log10(Te) against
 * log10(sigma*Te).
 *
 * @return a vector of Interp records covering all tabulated (Z, state)
 */
const vector<Interp> &interpolators() {
  static const vector<Interp> v = [] {
    vector<Interp> out;
    for (int Z = 1; Z <= 112; Z++) {
      const IchiharaEichlerData::Table *t =
          IchiharaEichlerData::FindTable(Z);
      for (int s = 0; s < (int)(sizeof(kStates) / sizeof(kStates[0])); s++) {
        Interp it;
        it.Z = Z;
        it.state = s;
        for (int i = 0; i < t->nRows; i++) {
          double sigma_barn = t->sigma[s * t->nRows + i];
          if (sigma_barn < 0.0)
            break; // column stops here for this state
          it.lx.push_back(log10(t->Te[i]));
          it.ly.push_back(log10(sigma_barn * t->Te[i]));
        }
        it.y2.resize(it.lx.size());
        spline((int)it.lx.size(), it.lx, it.ly, it.y2);
        out.push_back(move(it));
      }
    }
    return out;
  }();
  return v;
}

} // namespace

//////////////////////////////////////////////////////////////////////////
/**
 * @brief radiative-recombination cross section from the tabulated data of
 *        Ichihara & Eichler
 *
 * Natural cubic spline interpolation is performed in log10(E) against
 * log10(sigma*E).  Outside the range tabulated for the requested state the
 * cross section is returned as 0 rather than extrapolated.
 *
 * @param E_el_eV free-electron kinetic energy (eV)
 * @param Z       nuclear charge (1..112)
 * @param n       principal quantum number of the captured state
 * @param l       orbital angular momentum of the captured state
 * @param j       total angular momentum (half-integer, e.g. 1.5)
 *
 * @return sigma_RR (cm^2), or 0 if (Z, n, l, j) is not tabulated or
 *         E_el_eV falls outside the energy range for which the paper
 *         prints digits for that state
 */
double sigma_IchiharaEichlerRR(double E_el_eV, int Z, int n, int l, double j) {
  if (E_el_eV <= 0.0 || Z < 1 || Z > 112)
    return 0.0;

  int kappa = kappa_from_lj(l, j);
  int state = -1;
  for (int s = 0; s < (int)(sizeof(kStates) / sizeof(kStates[0])); s++) {
    if (kStates[s].n == n && kStates[s].kappa == kappa) {
      state = s;
      break;
    }
  }
  if (state < 0)
    return 0.0;

  for (auto &it : interpolators()) {
    if (it.Z != Z || it.state != state)
      continue;
    int npts = (int)it.lx.size();
    if (npts < 2)
      return 0.0;

    double lt = log10(E_el_eV);
    if (lt < it.lx.front() || lt > it.lx.back())
      return 0.0; // outside the range tabulated for this state
    double lsigmaE = splint(lt, npts, it.lx, it.ly, it.y2);
    double sigma_barn = pow(10.0, lsigmaE) / E_el_eV;
    return sigma_barn * 1e-24; // barn -> cm^2
  }
  return 0.0;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/**
 * @brief advance the radial Dirac continuum solution (P, Q) from r_a to r_b
 *
 * Fixed (kappa_c, c, E_el, z), r_b > r_a > 0.  Uses an adaptive-step,
 * error-controlled integrator (boost::numeric::odeint) when available;
 * otherwise falls back to a fine fixed-step RK4 sized to resolve both the
 * near-origin 1/r curvature and the asymptotic oscillation period 2*pi/k.
 *
 * @param kappa_c angular quantum number of the continuum state
 * @param c       speed of light in atomic units (1/alpha)
 * @param E_el    electron kinetic energy in Hartree
 * @param z       nuclear charge
 * @param r_a     starting radius
 * @param r_b     final radius
 * @param P       (in/out) large component, updated to r_b
 * @param Q       (in/out) small component, updated to r_b
 */
static void integrate_PQ(int kappa_c, double gamma_c, double c, double E_el,
                         double z, double r_a, double r_b, double &P,
                         double &Q) {
  if (r_b == r_a) return;

  // The raw radial functions P, Q span ~40 orders of magnitude between the
  // origin (P ~ r^gamma_c, down to ~1e-36 for gamma_c ~ 6, r0 = 1e-6) and
  // the far field (P ~ asymptotic amplitude, up to ~1e11), so no single
  // fixed absolute-error floor for the ODE stepper is small enough to
  // resolve the near-origin transient yet large enough to stay meaningful
  // at the far end.  Integrate whichever representation keeps the state
  // O(1)-ish at the current radius:
  //   * near the origin, where P, Q are tiny, integrate the rescaled
  //     variables p = P/r^gamma_c, q = Q/r^gamma_c (O(1) at r0 for every
  //     |kappa_c|), whose ODE reads
  //       dp/dr = -((κ+g)/r) p + (2c² + E + z/r) q/c
  //       dq/dr =  ((κ-g)/r) q - (E + z/r) p/c
  //   * once P, Q have grown well past the absolute floor (so the relative
  //     tolerance governs), continue in the raw variables, which stay
  //     O(amplitude) in the far field instead of collapsing like 1/r^gamma.
  // This avoids both the under-resolved near-origin transient (which
  // corrupts the far-field phase/amplitude and showed up as runaway DIR
  // cross sections for l >= 4) and the under-resolved far field (the
  // rescaled state falls back below the absolute floor as p ~ A/r^gamma).
  const double thr = 1e-25; // ~1e5 above the abs_err=1e-30 floor
  const bool use_raw = (fabs(P) > thr) || (fabs(Q) > thr);

  double rp_a = use_raw ? 1.0 : pow(r_a, gamma_c);
  std::array<double, 2> y{use_raw ? P : P / rp_a,
                          use_raw ? Q : Q / rp_a};

  auto deriv = [&](const std::array<double, 2> &yv, std::array<double, 2> &dy,
                    double rr) {
    if (use_raw) {
      dy[0] = -(kappa_c / rr) * yv[0] + (2.0 * c * c + E_el + z / rr) * yv[1] / c;
      dy[1] =  (kappa_c / rr) * yv[1] - (E_el + z / rr) * yv[0] / c;
    } else {
      dy[0] = -((kappa_c + gamma_c) / rr) * yv[0] +
              (2.0 * c * c + E_el + z / rr) * yv[1] / c;
      dy[1] =  ((kappa_c - gamma_c) / rr) * yv[1] -
              (E_el + z / rr) * yv[0] / c;
    }
  };

#ifdef HYDROCAL_HAVE_BOOST_ODEINT
  namespace odeint = boost::numeric::odeint;
  double abs_err = 1e-30, rel_err = 1e-10;
  auto stepper = odeint::make_controlled(
      abs_err, rel_err, odeint::runge_kutta_dopri5<std::array<double, 2>>());
  // The initial step must also respect the local curvature scale ~r_a, not
  // just the requested span: for the very first marching call out of r0
  // (~1e-6), a span/100 guess is many orders of magnitude larger than r0
  // itself, so the stepper's first step lands well past the near-origin
  // region (steepest for high |kappa_c|, gamma_c ~ 2+) before the error
  // controller has adapted -- producing a silently-inaccurate first step
  // whose error looks small only because it's being compared against a
  // relative tolerance on values that are themselves still tiny there.
  double h0 = std::min((r_b - r_a) / 100.0, std::max(r_a, 1e-12) * 0.05);
  odeint::integrate_adaptive(stepper, deriv, y, r_a, r_b, h0);
#else
  double k = sqrt(E_el * (E_el + 2.0 * c * c)) / c;
  double span = r_b - r_a;
  double steps_for_oscillation = 50.0 * span * k;
  double steps_for_curvature = 50.0 * span / std::max(r_a, 1e-9);
  int nsub = (int)std::max({20.0, steps_for_oscillation, steps_for_curvature});
  nsub = std::min(nsub, 500000);
  double h = span / nsub;
  double rr = r_a;
  for (int i = 0; i < nsub; i++) {
    std::array<double, 2> k1, k2, k3, k4, ytmp;
    deriv(y, k1, rr);
    ytmp = {y[0] + 0.5 * h * k1[0], y[1] + 0.5 * h * k1[1]};
    deriv(ytmp, k2, rr + 0.5 * h);
    ytmp = {y[0] + 0.5 * h * k2[0], y[1] + 0.5 * h * k2[1]};
    deriv(ytmp, k3, rr + 0.5 * h);
    ytmp = {y[0] + h * k3[0], y[1] + h * k3[1]};
    deriv(ytmp, k4, rr + h);
    y[0] += h / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
    y[1] += h / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
    rr += h;
  }
#endif

  double rp_b = use_raw ? 1.0 : pow(r_b, gamma_c);
  P = rp_b * y[0];
  Q = rp_b * y[1];
}



//////////////////////////////////////////////////////////////////////////
/**
 * @brief Parameters of the closed-form continuum construction
 *
 * used in dirac_continuum_PQ_closed
 *
 * Everything here depends only on (E_el_eV, z, kappa_c), so it
 * is cached across the per-grid-point calls, and it is shared verbatim by
 * the analytic bound-free integral (dirac_rr_bf_analytic), guaranteeing
 * that the analytic and numerical paths use the same continuum
 */
struct ContinuumParams {
  double k;   ///< wavenumber of the electron 
  double s;   ///< sqrt(kappa_c² - Z² alpha²);
  double b1;  ///< (E_el + q0 * q2 * W0) / (c * (q2 - q0));
  double c12; ///< (E_el + q2 * q2 * W0) / (c * (q2 - q0));
  double q0;  ///< (s + kappa_c) * c / Z;
  double q2;  ///< (kappa_c - s) * c / Z;
  double N;   ///< energy normalisation, fixed by the exact asymptotic branch amplitudes.
  complex<double> a; ///< s*[1 + i*(b1/k)];
};

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Initialization of the parameters of the closed-form continuum construction
 *
 * @param E_el_eV electron energy in eV
 * @param z nuclear charge
 * @param kappa_c relativistic angular momentum quantum number of the continuum electron
 */
static ContinuumParams dirac_continuum_params(double E_el_eV, double z, int kappa_c) {
  using hydroconst::alpha;
  using hydroconst::pi;
  static thread_local double cache_E = numeric_limits<double>::quiet_NaN();
  static thread_local double cache_z = numeric_limits<double>::quiet_NaN();
  static thread_local int cache_kappa = 0;
  static thread_local ContinuumParams cp;
  bool hit = (cache_E == E_el_eV && cache_z == z && cache_kappa == kappa_c);
  if (!hit) {
    double c = 1.0 / alpha;
    double E_el = E_el_eV / (2.0 * hydroconst::Ryd_eV); // Hartree
    cp.k = sqrt(E_el * (E_el + 2.0 * c * c)) / c;
    double za = z * alpha;
    cp.s = sqrt(kappa_c * kappa_c - za * za);
    double W0 = 2.0 * c * c + E_el;
    cp.q0 = (cp.s + kappa_c) * c / z;
    cp.q2 = (kappa_c - cp.s) * c / z;
    cp.b1 = (E_el + cp.q0 * cp.q2 * W0) / (c * (cp.q2 - cp.q0));
    cp.c12 = (E_el + cp.q2 * cp.q2 * W0) / (c * (cp.q2 - cp.q0));
    cp.a = complex<double>(cp.s, cp.s * cp.b1 / cp.k);
    double A_asymp = sqrt(2.0 / (pi * cp.k)) * sqrt(1.0 + E_el / (c * c));

    // Exact asymptotic branch amplitudes of u1(r) ~ C_in r^{a-2s} e^{-ikr}
    // + C_out r^{-a} e^{ikr}, and the branch-ratio connection coefficients
    // t_in, t_out that relate u2 to u1 in the far field.  |C_in|, |C_out|
    // are O(1) even at low energy (the e^{-pi|y|/2} of |Gamma(a)| with
    // y = Im(a) ~ -pi*eta is cancelled by the e^{+pi|y|/2} of |(-2ik)^{a-2s}|),
    // but each factor underflows separately, so the product is formed as
    // exp(sum of logs) via clngamma_cmplx.
    complex<double> twos = 2.0 * cp.s;
    complex<double> log_minus2ik = log(complex<double>(0.0, -2.0 * cp.k));
    complex<double> log_2ik = log(complex<double>(0.0, 2.0 * cp.k));
    complex<double> C_in =
        cgamma_cmplx(twos) *
        exp(-clngamma_cmplx(cp.a) + (cp.a - twos) * log_minus2ik);
    complex<double> C_out =
        cgamma_cmplx(twos) *
        exp(-clngamma_cmplx(twos - cp.a) - cp.a * log_2ik);
    complex<double> t_in = -(complex<double>(0.0, cp.k) + cp.b1) / cp.c12;
    complex<double> t_out = (complex<double>(0.0, cp.k) - cp.b1) / cp.c12;
    complex<double> CinP = (1.0 + t_in) * C_in;
    complex<double> CoutP = (1.0 + t_out) * C_out;
    complex<double> CinQ = (cp.q0 + cp.q2 * t_in) * C_in;
    complex<double> CoutQ = (cp.q0 + cp.q2 * t_out) * C_out;
    double w = norm(CinP) + norm(CoutP) + norm(CinQ) + norm(CoutQ);
    cp.N = (w > 1e-300) ? A_asymp / sqrt(w) : 0.0;
    cache_E = E_el_eV;
    cache_z = z;
    cache_kappa = kappa_c;
  }
  return cp;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief energy-normalised Dirac continuum wavefunction P(r), Q(r) in closed form
 *
 * The continuum wavefunction of the radial Dirac equation is written in
 * closed form in terms of the confluent hypergeometric function 1F1
 * (sec:closedform of DiracTransitionRates.tex):
 *   P_c(r) = N r^s Re[u1(r) + u2(r)]
 *   Q_c(r) = N r^s Re[q0 u1(r) + q2 u2(r)]
 * with u1 = e^{ikr} 1F1(a; 2s; -2ikr), a = s(1 + i b1/k), and
 * u2 = (u1' - b1 u1)/c12.  The normalisation N is fixed by the exact
 * asymptotic branch amplitudes C_in, C_out of u1 (via the connection
 * coefficients t_in, t_out), so that the far field carries the
 * energy-normalised amplitude A_asymp = sqrt(2/(pi k)) sqrt(1 + E_el/c^2).
 * Unlike dirac_continuum_PQ this is an O(1) point evaluation that needs no
 * integration and suffers no deep-cancellation loss when the 1F1 is
 * evaluated with the exact engine.
 *
 * @param E_el_eV electron kinetic energy in eV
 * @param z       nuclear charge
 * @param kappa_c angular quantum number of the continuum state
 * @param r       radius (a.u.) at which the wavefunction is requested
 * @param P       (output) large component at radius r
 * @param Q       (output) small component at radius r
 */
static void dirac_continuum_PQ_closed(double E_el_eV, double z, int kappa_c,
                                      double r, double &P, double &Q) {
  ContinuumParams cp = dirac_continuum_params(E_el_eV, z, kappa_c);
  if (cp.k < 1e-30) {
    P = 0.0;
    Q = 0.0;
    return;
  }

  double kr = cp.k * r;
  complex<double> eikr(cos(kr), sin(kr));
  complex<double> f, g;
  hypconfl_cmplx(cp.a.real(), cp.a.imag(), 2.0 * cp.s, 0.0, -2.0 * kr, f);
  hypconfl_cmplx(cp.a.real() + 1.0, cp.a.imag(), 2.0 * cp.s + 1.0, 0.0,
                 -2.0 * kr, g);
  // u1' = ik e^{ikr} [1F1(a;2s;z) - (a/s) 1F1(a+1;2s+1;z)]
  complex<double> u1 = eikr * f;
  complex<double> u1p = complex<double>(0.0, cp.k) * eikr *
                        (f - (cp.a / cp.s) * g);
  complex<double> u2 = (u1p - cp.b1 * u1) / cp.c12;

  double rpow = pow(r, cp.s) * cp.N;
  P = rpow * real(u1 + u2);
  Q = rpow * real(cp.q0 * u1 + cp.q2 * u2);
}

////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Parameters of the hydrogenic bound state (n, kappa)
 *
 * Shared by the analytic and numerical bound-free integrals.
 *
 * The a_coef, b_coef already carry the (2 lam)^i scale.
 */
struct BoundWaveParams {
  double gamma; ///< first amplitude factor
  double eps; ///< second amplitude factor
  double lam; ///< decay constant
  double two_lam;
  double C;  ///< nomalization constant
  int nr;
  vector<double> a_coef;  ///< coefficients in P_b = C sqrt(1+eps) (2 lam r)^gamma e^{-lam r} sum_i a_coef[i] r^i
  vector<double> b_coef;  ///< coefficients in Q_b = C sqrt(1-eps) (2 lam r)^gamma e^{-lam r} sum_i b_coef[i] r^i
};


////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Initialization of the parameters of the hydrogenic bound state (n, kappa)
 *
 * @param z nuclear charge
 * @param n principal qauntum number
 * @param kappa relativistic angular momentum quantum number 
 */
 static BoundWaveParams dirac_bound_wave_params(double z, int n, int kappa) {
  using hydroconst::alpha;
  BoundWaveParams bw;
  double za = z * alpha;
  bw.nr = n - abs(kappa);
  bw.gamma = sqrt(kappa * kappa - za * za);
  double Nprime = bw.nr + bw.gamma;
  bw.eps = 1.0 / sqrt(1.0 + za * za / (Nprime * Nprime));
  bw.lam = z / sqrt(Nprime * Nprime + za * za);   // a.u.
  bw.two_lam = 2.0 * bw.lam;
  bw.C = dirac_normalisation(z, n, kappa);
  build_FG_coeffs(kappa, bw.nr, bw.gamma, bw.eps, bw.two_lam, z,
                  bw.a_coef, bw.b_coef);
  return bw;
}

///////////////////////////////////////////////////////////////////////////////////
/**
 * @brief  Analytic bound-free radial integral for E1 transition
 *
 *  for the bound state (n, kappa) and the E1-coupled continuum channel
 *  kappa_c, evaluated without numerical integration.  With the closed-form
 * continuum
 *   P_c = N r^s Re[u1+u2],   Q_c = N r^s Re[q0 u1 + q2 u2],
 *   u1 = e^{ikr} 1F1(a;2s;-2ikr),   u2 = (u1' - b1 u1)/c12,
 * every polynomial term of the integrand is of the form
 *   r^{gamma+s+1+i} e^{-(lam-ik)r} 1F1(a'; b'; -2ikr),
 * which integrates term-by-term via the Laplace transform of 1F1
 * (DLMF 13.10.7),
 *   int_0^inf t^{nu-1} e^{-mu t} 1F1(a';b';kap t) dt
 *       = Gamma(nu) mu^{-nu} 2F1(a', nu; b'; kap/mu),   Re nu > 0,
 * with mu = lam - ik and kap/mu = -2ik/(lam-ik).  R_bf is therefore a
 * finite sum of Gamma(nu) mu^{-nu} 2F1(...) values -- the relativistic
 * analog of the (nonrelativistic) Gordon element: exact at every energy
 * and free of the deep-cancellation loss that limits the numerical
 * quadrature at high E_kin.  Returns false if a hypergeometric evaluation
 * was not trustworthy (portable 2F1 only), in which case the caller falls
 * back to the numerical path.
 *
 * @param E_el_eV electron kinetic energy in eV
 * @param z       nuclear charge
 * @param kappa_c angular quantum number of the continuum state
 * @param bw parameters describing the bound-state wave function
 * @param on exit: R_bf = int_0^inf r [P_b P_c + Q_b Q_c] dr
 *
 * @return true/false if evaluation was successful/unsuccessful 
 */
static bool dirac_rr_bf_analytic(double E_el_eV, double z, int kappa_c,
                                 const BoundWaveParams &bw, double &R_bf) {
  ContinuumParams cp = dirac_continuum_params(E_el_eV, z, kappa_c);
  if (cp.k < 1e-30) {
    R_bf = 0.0;
    return true;
  }

  complex<double> ik(0.0, cp.k);
  complex<double> mu(bw.lam, -cp.k);        // lam - i k
  complex<double> z2f1 = (-2.0 * ik) / mu;  // -2ik/(lam-ik)
  complex<double> logmu = log(mu);

  // Branch coefficients of u1+u2 and q0 u1 + q2 u2 in f = 1F1(a;2s;z) and
  // g = 1F1(a+1;2s+1;z), with u2 = (u1' - b1 u1)/c12 and
  // u1' = ik e^{ikr}[f - (a/s) g]:
  //   u1+u2      = e^{ikr}[ (1 + (ik-b1)/c12) f - (ik a/s)/c12 g ]
  //   q0 u1+q2 u2= e^{ikr}[ (q0 + q2(ik-b1)/c12) f - (q2 ik a/s)/c12 g ]
  complex<double> c_f0 = (ik - cp.b1) / cp.c12;
  complex<double> c_g0 = -(ik * cp.a / cp.s) / cp.c12;
  complex<double> cQ_f = cp.q0 + cp.q2 * c_f0;
  complex<double> cQ_g = cp.q2 * c_g0;

  complex<double> total(0.0, 0.0);
  for (int i = 0; i <= bw.nr; i++) {
    double nu = bw.gamma + cp.s + i + 2;
    complex<double> gpref = exp(std::lgamma(nu) - nu * logmu);
    complex<double> F1, F2;
    if (!hypgeo2_cmplx(cp.a, complex<double>(nu, 0.0),
                       complex<double>(2.0 * cp.s, 0.0), z2f1, F1))
      return false;
    if (!hypgeo2_cmplx(cp.a + 1.0, complex<double>(nu, 0.0),
                       complex<double>(2.0 * cp.s + 1.0, 0.0), z2f1, F2))
      return false;
    complex<double> ap = sqrt(1.0 + bw.eps) * bw.a_coef[i];
    complex<double> aq = sqrt(1.0 - bw.eps) * bw.b_coef[i];
    complex<double> Aamp = bw.C * (ap * (1.0 + c_f0) + aq * cQ_f);
    complex<double> Bamp = bw.C * (ap * c_g0 + aq * cQ_g);
    total += (Aamp * F1 + Bamp * F2) * gpref;
  }
  R_bf = cp.N * exp(bw.gamma * log(bw.two_lam)) * real(total);
  return true;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief energy-normalised Dirac continuum wavefunction P(r), Q(r)
 *
 * The radial Dirac equation is integrated numerically for the given
 * kinetic energy E_el_eV and quantum number kappa_c.  The solution is
 * normalised so that
 *   int_0^inf [P_E P_{E'} + Q_E Q_{E'}] dr = delta(E - E').
 *
 * @param E_el_eV electron kinetic energy in eV
 * @param z       nuclear charge
 * @param kappa_c angular quantum number of the continuum state
 * @param r       radius (a.u.) at which the wavefunction is requested
 * @param P       (output) large component at radius r
 * @param Q       (output) small component at radius r
 */
static void dirac_continuum_PQ(double E_el_eV, double z, int kappa_c,
                                double r, double &P, double &Q) {
  using hydroconst::alpha;

  double c = 1.0 / alpha;
  double E_el = E_el_eV / (2.0 * hydroconst::Ryd_eV); // Hartree

  // Relativistic wave number
  double k = sqrt(E_el * (E_el + 2.0 * c * c)) / c;
  if (k < 1e-30) { P = 0.0; Q = 0.0; return; }

  double za = z * alpha;
  double gamma_c = sqrt(kappa_c * kappa_c - za * za);

  // Integrate the radial Dirac ODE from r0 outward to r (if r > r0).
  // ODE:  dP/dr = -(κ/r) P + (E + c² + Z/r) Q / c
  //       dQ/dr =  (κ/r) Q - (E - c² + Z/r) P / c
  //
  // For a continuum state: E = c² + E_el
  //   dP/dr = -(κ/r) P + (2c² + E_el + Z/r) Q / c
  //   dQ/dr =  (κ/r) Q - (E_el + Z/r) P / c

  double r0 = 1e-6;

  // Starting values at r0 from a two-term Frobenius series of the coupled
  // ODE above (P = r^s (1 + p1 r), Q = r^s (q0 + q1 r), s = gamma_c). The
  // leading-order ratio q0 = c(s+kappa_c)/z (sign convention cross-checked
  // against the bound-state series in build_FG_coeffs: bu[0]=gamma+kappa)
  // degenerates toward 0 as z*alpha -> 0 for every kappa_c < 0 (since then
  // gamma_c -> |kappa_c| = -kappa_c), which makes the leading-order-only
  // starting condition numerically invalid for those channels at low Z.
  // The next-order term p1, q1 (from matching the O(r^s) terms in the ODE)
  // is well-conditioned in that limit (det = 2*gamma_c+1 != 0 always) and
  // removes the degeneracy.
  double s = gamma_c;
  double W = 2.0 * c * c + E_el;
  double q0 = (s + kappa_c) * c / z;
  double Acoef = s + 1.0 + kappa_c;
  double Ccoef = s + 1.0 - kappa_c;
  double Bcoef = z / c;
  double det = 2.0 * s + 1.0;
  double p1 = (W * q0 * Ccoef - Bcoef * E_el) / (c * det);
  double q1 = -(Acoef * E_el + Bcoef * W * q0) / (c * det);

  double rpow = pow(r0, gamma_c);
  double P0 = rpow * (1.0 + p1 * r0);
  double Q0 = rpow * (q0 + q1 * r0);

  // Normalise to the energy-normalised asymptotic amplitude
  // At large r: P ~ A sin(kr - lπ/2 + δ + η ln(2kr))
  // Energy normalisation gives: A = sqrt(2/(πk)) * sqrt(1 + E_el/c²)
  double A_asymp = sqrt(2.0 / (hydroconst::pi * k)) * sqrt(1.0 + E_el / (c * c));

  if (r <= r0) { P = P0 * A_asymp; Q = Q0 * A_asymp; return; }

  // The normalisation constant only depends on (E_el_eV, z, kappa_c), not
  // on r, but finding it requires integrating well past r into the
  // asymptotic region. dirac_rr_xsec() calls this function many times
  // (once per radial-grid point) with the same (E_el_eV, z, kappa_c), so
  // cache the expensive asymptotic search and reuse it for every r.
  static thread_local double cached_E = std::numeric_limits<double>::quiet_NaN();
  static thread_local double cached_z = std::numeric_limits<double>::quiet_NaN();
  static thread_local int cached_kappa = 0;
  static thread_local double cached_norm = 1.0;
  // Marching state: the raw (unnormalised) solution at the last radius this
  // function was asked to evaluate for the same (kappa_c, z, E_el_eV).  RR
  // cross-section calculations sweep r over a large radial grid in ascending
  // order, and each grid step is tiny compared to the ~500 oscillation
  // periods separating r0 from the outer grid edge at high energy;
  // integrating from r0 for every single grid point forces the adaptive
  // stepper to cross that whole span in one call (with only a crude initial
  // step-size guess) each time, which both wastes time and, empirically,
  // degrades accuracy enough to corrupt the bound-free radial integral at
  // high energy. Continuing from the last visited radius instead keeps each
  // integrate_PQ() call short and well resolved.
  static thread_local double cached_r_last = std::numeric_limits<double>::quiet_NaN();
  static thread_local double cached_P_last = 0.0;
  static thread_local double cached_Q_last = 0.0;

  bool cache_hit =
      (cached_kappa == kappa_c && cached_z == z && cached_E == E_el_eV);

  if (!cache_hit) {
    // Integrate outward period by period, tracking the RMS envelope
    // amplitude sqrt(P^2+Q^2) sampled across each period, until it stops
    // changing (envelope has reached its asymptotic, slowly-varying
    // value) or a generous safety cap on the number of periods is hit.
    //
    // The nominal oscillation period 2pi/k is only the r->infinity limit.
    // At finite r the attractive Coulomb tail z/r adds a slowly-varying
    // contribution to the local phase advance, dtheta/dr = k + eta/r (eta
    // is the relativistic Sommerfeld parameter below, matching the
    // 1+E_el/c^2 factor already in A_asymp: the electron speeds up as it
    // falls toward the nucleus, so the local wavenumber is *higher* than
    // k, not lower). Sampling nsamp_per_period points across a *fixed*
    // period 2pi/k therefore drifts out of phase with the true oscillation
    // by O(eta/(k r)) per period; that mismatch is a systematic (not
    // random) bias in the RMS envelope that only falls off as 1/r, so the
    // old fixed-period sampling needed impractically many periods to reach
    // even 1% accuracy for large eta = Z*(1+E_el/c^2)/k (high Z, low
    // energy). Using the local period 2pi/(k+eta/r) instead removes this
    // bias at its source.
    double eta = z * (1.0 + E_el / (c * c)) / k;

    // Below r ~ eta/k the 1/r correction term would formally exceed k
    // itself; this is exactly the near-origin region where the solution
    // isn't a simple k-oscillation yet anyway (it is still relaxing out of
    // the r^gamma_c Frobenius regime), so cap the correction there at its
    // r=eta/k value instead of letting it diverge.
    double r_floor = eta / k;

    // Even with the corrected local period, the envelope still shows a
    // residual slow beat (neighbouring periods can coincidentally agree to
    // within a tight tolerance while the underlying trend is still
    // drifting), so comparing single consecutive periods is not a
    // reliable convergence test. Instead average over a window of periods
    // and require the windowed average itself to stop moving between one
    // window and the next -- this is insensitive to the beat and to
    // sample-to-sample noise.
    double rr = r0, Pw = P0, Qw = Q0;
    double amp_rms = 0.0;
    const int nsamp_per_period = 8;
    const int max_periods = 20000;
    // tol=5e-4, window=200 was chosen empirically: at the worst case tested
    // (Z=112, kappa_c=1, E_el=10 eV) it reaches ~90% of the way from the
    // original tol=5e-3/window=100 result to the value obtained by an
    // (expensive, 50x more periods) tol=1e-5/window=1000 cross-check, while
    // still converging in a few thousand periods rather than tens of
    // thousands. A genuine residual (~1 percentage point at that worst
    // case) remains uncaptured even by the tight cross-check -- see
    // sec:ichiharavalidation in DiracTransitionRates.tex for the
    // Z-dependent trend this leaves behind and why it is not believed to be
    // a further convergence artifact.
    const double tol = 5e-4;
    // A single passing window-vs-window comparison is still not reliable:
    // during the still-rising, noisy climb toward the true plateau the
    // windowed mean can pass through a coincidental flat spot (two
    // adjacent windows agreeing to within tol while the mean a few hundred
    // periods later is still substantially higher) -- observed directly
    // for high Z at low energy. Guard against this by requiring several
    // *consecutive* passing comparisons (i.e. the mean must stay flat for
    // consec_needed*window periods in a row), not just one.
    const int window = 200;
    const int consec_needed = 5;
    int consec = 0;
    std::vector<double> hist;
    hist.reserve(max_periods);

    for (int iper = 0; iper < max_periods; iper++) {
      double k_local = k + eta / std::max(rr, r_floor);
      double period = 2.0 * hydroconst::pi / k_local;
      double sumsq = 0.0;
      for (int is = 0; is < nsamp_per_period; is++) {
        double r_next = rr + period / nsamp_per_period;
        integrate_PQ(kappa_c, gamma_c, c, E_el, z, rr, r_next, Pw, Qw);
        rr = r_next;
        sumsq += Pw * Pw + Qw * Qw;
      }
      amp_rms = sqrt(sumsq / nsamp_per_period);
      hist.push_back(amp_rms);
      int n = (int)hist.size();
      if (n >= 2 * window) {
        double sum_prev = 0.0, sum_cur = 0.0;
        for (int j = n - 2 * window; j < n - window; j++) sum_prev += hist[j];
        for (int j = n - window; j < n; j++) sum_cur += hist[j];
        double mean_prev = sum_prev / window;
        double mean_cur = sum_cur / window;
        if (fabs(mean_cur - mean_prev) < tol * mean_cur) {
          consec++;
          amp_rms = mean_cur;
          if (consec >= consec_needed) break;
        } else {
          consec = 0;
        }
      }
    }

    cached_norm = (amp_rms > 1e-300) ? A_asymp / amp_rms : 1.0;
    cached_E = E_el_eV;
    cached_z = z;
    cached_kappa = kappa_c;
    cached_r_last = r0;
    cached_P_last = P0;
    cached_Q_last = Q0;
  }

  double r_from = r0, Pr = P0, Qr = Q0;
  if (cache_hit && r >= cached_r_last) {
    r_from = cached_r_last;
    Pr = cached_P_last;
    Qr = cached_Q_last;
  }
  integrate_PQ(kappa_c, gamma_c, c, E_el, z, r_from, r, Pr, Qr);
  cached_r_last = r;
  cached_P_last = Pr;
  cached_Q_last = Qr;

  P = Pr * cached_norm;
  Q = Qr * cached_norm;
}

////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Fully-retarded (plane-wave) bound-free matrix elements
 *
 * The retarded operator e^{-ik_ph r cosT} replaces the dipole (long-wavelength) r factor:
 *  expanding the plane wave in spherical Bessel functions,
 *   e^{-ik_ph r cosT} = sum_L (2L+1) (-i)^L j_L(k_ph r) P_L(cosT),
 * gives the radial integrals
 *  - Rp[L] = int_0^inf P_c(r) Q_b(r) j_L(k_ph r) dr,
 *  - Rm[L] = int_0^inf Q_c(r) P_b(r) j_L(k_ph r) dr,
 * summed over all L with the angular coefficients of ang_int.  Both
 * parities and all multipoles are allowed (the angular coefficients
 * decide), so the channel set |j_c - j_b| = 0, 1, 2, ... is complete.
 *
 *Analytic (retarded) bound-free radial integrals Rp[L], Rm[L] in closed
 * form, again via the Laplace transform of 1F1 (DLMF 13.10.7).  With
 *   j_L(kp r) = (1/2) sum_{m=0}^L C_{Lm} (kp r)^(-m-1)
 *                   [ i^{m-L-1} e^{+ikp r} + (-i)^{m-L-1} e^{-ikp r} ],
 *   C_{Lm} = (L+m)! / (2^m m! (L-m)!),
 * every polynomial term of P_c Q_b (resp. Q_c P_b) is
 *   r^{gamma+s+i-m-1} e^{-mu r} 1F1(a'; b'; -2ikc r)
 * with nu = gamma + s + i - m and the two branches
 *   mu_plus  = lam - i(kc + kp)   (the e^{+ikp r} branch),
 *   mu_minus = lam - i(kc - kp)   (the e^{-ikp r} branch),
 * integrating to Gamma(nu) mu^{-nu} 2F1(a', nu; b'; -2ikc/mu).  Rp[L] and
 * Rm[L] are therefore finite sums of Gamma(nu) mu^{-nu} 2F1(...) values,
 * the retarded generalisation of dirac_rr_bf_analytic: exact at every
 * energy, all multipoles, no quadrature.  Returns false if a
 * hypergeometric evaluation was not trustworthy (portable 2F1 only), in
 * The retarded radial integrals
 *   Rp[L] = int_0^inf P_c Q_b j_L(kp r) dr,  Rm[L] = int_0^inf Q_c P_b j_L(kp r) dr
 * in closed form via the branch decomposition of j_L(kp r).  For photon
 * wavenumbers kp below ~0.5 the m-branch terms of order kp^{-(m+1)} cancel
 * to leave j_L ~ (kp r)^L at a depth that exceeds double precision, and
 * Gamma(nu) approaches a pole as nu = gamma + s + i - m -> 0 (m ~ s at low Z).
 * Such L values are flagged in 'unreliable' and the caller replaces
 * them with the small-kp Taylor series (dirac_rr_bf_retarded_taylor), which
 *case the caller falls back to dirac_rr_bf_retarded_numeric.
 */
static bool dirac_rr_bf_retarded(double E_el_eV, double z, int kappa_c,
                                 const BoundWaveParams &bw, double kp,
                                 int Lmax, vector<double> &Rp,
                                 vector<double> &Rm,
                                 vector<char> &unreliable) {
  ContinuumParams cp = dirac_continuum_params(E_el_eV, z, kappa_c);
  unreliable.assign(Lmax + 1, 0);
  if (cp.k < 1e-30) {
    Rp.assign(Lmax + 1, 0.0);
    Rm.assign(Lmax + 1, 0.0);
    return true;
  }

  complex<double> ik(0.0, cp.k);
  complex<double> mu_plus(bw.lam, -(cp.k + kp));   // lam - i(kc+kp)
  complex<double> mu_minus(bw.lam, -(cp.k - kp));  // lam - i(kc-kp)
  complex<double> z2_plus = (-2.0 * ik) / mu_plus;
  complex<double> z2_minus = (-2.0 * ik) / mu_minus;
  complex<double> logmu_plus = log(mu_plus);
  complex<double> logmu_minus = log(mu_minus);

  // Branch coefficients of u1+u2 and q0 u1 + q2 u2 in f = 1F1(a;2s;z) and
  // g = 1F1(a+1;2s+1;z), as in dirac_rr_bf_analytic.
  complex<double> c_f0 = (ik - cp.b1) / cp.c12;
  complex<double> c_g0 = -(ik * cp.a / cp.s) / cp.c12;
  complex<double> cQ_f = cp.q0 + cp.q2 * c_f0;
  complex<double> cQ_g = cp.q2 * c_g0;

  // (1i)^n for the j_L branch phases; the e^{-ikp r} branch uses
  // (-1i)^n = conj((1i)^n).  Period 4, valid for negative n as well.
  auto ipow = [](int n) {
    int r = n % 4;
    if (r < 0) r += 4;
    static const complex<double> t[4] = {complex<double>(1.0, 0.0),
                                         complex<double>(0.0, 1.0),
                                         complex<double>(-1.0, 0.0),
                                         complex<double>(0.0, -1.0)};
    return t[r];
  };

  // Log-factorials up to 2*Lmax for the C_{Lm} coefficients.
  vector<double> lfact(2 * Lmax + 1);
  lfact[0] = 0.0;
  for (size_t k = 1; k < lfact.size(); k++)
    lfact[k] = lfact[k - 1] + log((double)k);
  double logkp = log(kp);

  // The hypergeometric values depend only on (i, m, branch), not on L (nu
  // is shared by every L >= m): cache them once per (i, m) and reuse.
  struct HTerm {
    complex<double> gp, F1, F2;
  };
  int ntab = (bw.nr + 1) * (Lmax + 1);
  vector<HTerm> term_plus(ntab), term_minus(ntab);
  for (int i = 0; i <= bw.nr; i++) {
    for (int m = 0; m <= Lmax; m++) {
      double nu = bw.gamma + cp.s + i - m;
      complex<double> cnu(nu, 0.0);
      int idx = i * (Lmax + 1) + m;
      // mpmath.loggamma(x) for x < 0 carries the branch phase
      // e^{-i pi floor(x)} that reproduces sign(Gamma(x)); std::lgamma
      // returns ln|Gamma(x)|, so fold the sign in explicitly.
      double lgnu = std::lgamma(nu);
      double gsign = (nu >= 0.0) ? 1.0 : (((int)floor(nu) % 2) ? -1.0 : 1.0);
      term_plus[idx].gp = gsign * exp(lgnu - nu * logmu_plus);
      term_minus[idx].gp = gsign * exp(lgnu - nu * logmu_minus);
      if (!hypgeo2_cmplx(cp.a, cnu, complex<double>(2.0 * cp.s, 0.0),
                         z2_plus, term_plus[idx].F1))
        return false;
      if (!hypgeo2_cmplx(cp.a + 1.0, cnu,
                         complex<double>(2.0 * cp.s + 1.0, 0.0), z2_plus,
                         term_plus[idx].F2))
        return false;
      if (!hypgeo2_cmplx(cp.a, cnu, complex<double>(2.0 * cp.s, 0.0),
                         z2_minus, term_minus[idx].F1))
        return false;
      if (!hypgeo2_cmplx(cp.a + 1.0, cnu,
                         complex<double>(2.0 * cp.s + 1.0, 0.0), z2_minus,
                         term_minus[idx].F2))
        return false;
    }
  }

  double pref = cp.N * bw.C * exp(bw.gamma * log(bw.two_lam));
  double s1p = pref * sqrt(1.0 - bw.eps);  // Rp scale (Q_b)
  double s1m = pref * sqrt(1.0 + bw.eps);  // Rm scale (P_b)
  // Bound-state density is negligible beyond r ~ 2 N'/lam (N' = nr + gamma
  // the effective principal quantum number): there the j_L(kp r) turning
  // point L/kp lies beyond the wavefunction support and the multipole
  // integral is exponentially small.  In that regime the branch m-sum can
  // still evaluate to a spurious moderate-depth value (depth ~ 1e6..1e13,
  // below the double-precision limit but orders above the true magnitude),
  // so those multipoles are flagged whenever the sum shows any residual
  // cancellation.
  double r_bound = 2.0 * (bw.nr + bw.gamma) / bw.lam;

  Rp.assign(Lmax + 1, 0.0);
  Rm.assign(Lmax + 1, 0.0);
  for (int L = 0; L <= Lmax; L++) {
    // Kahan-compensated (two-sum) accumulation of the alternating m-sum.
    complex<double> total_p(0.0, 0.0), total_m(0.0, 0.0);
    complex<double> cpk(0.0, 0.0), cmk(0.0, 0.0);
    double peak_p = 0.0, peak_m = 0.0;
    for (int i = 0; i <= bw.nr; i++) {
      for (int m = 0; m <= L; m++) {
        int expn = m - L - 1;
        double logC = lfact[L + m] - m * log(2.0) - lfact[m] - lfact[L - m];
        complex<double> half = 0.5 * exp(logC - (m + 1.0) * logkp);
        complex<double> coef_plus = half * ipow(expn);
        complex<double> coef_minus = half * conj(ipow(expn));
        int idx = i * (Lmax + 1) + m;
        const HTerm &tp = term_plus[idx];
        const HTerm &tm = term_minus[idx];
        complex<double> pterm =
            bw.b_coef[i] *
            (coef_plus * tp.gp * ((1.0 + c_f0) * tp.F1 + c_g0 * tp.F2) +
             coef_minus * tm.gp * ((1.0 + c_f0) * tm.F1 + c_g0 * tm.F2));
        complex<double> mterm =
            bw.a_coef[i] *
            (coef_plus * tp.gp * (cQ_f * tp.F1 + cQ_g * tp.F2) +
             coef_minus * tm.gp * (cQ_f * tm.F1 + cQ_g * tm.F2));
        peak_p = max(peak_p, abs(pterm));
        peak_m = max(peak_m, abs(mterm));
        complex<double> y = pterm - cpk;
        complex<double> t = total_p + y;
        cpk = (t - total_p) - y;
        total_p = t;
        y = mterm - cmk;
        t = total_m + y;
        cmk = (t - total_m) - y;
        total_m = t;
      }
    }
    double RpL = s1p * real(total_p);
    double RmL = s1m * real(total_m);
    // The m-branch expansion represents j_L(kp r) as a difference of terms of
    // order kp^{-(m+1)}; the cancellation depth peak/|raw sum| must stay far
    // below double precision.  (The comparison is against the un-prefactored
    // sums -- s1p/s1m scale the result by up to ~1e4, which would otherwise
    // hide the cancellation.)  Multipoles beyond the turning-point boundary
    // are additionally screened at a lower depth to catch the spurious
    // moderate-cancellation values they would otherwise inject.  Flagged
    // multipoles are patched with the Taylor series by the caller.
    double depth_p = peak_p / (abs(total_p) + 1e-300);
    double depth_m = peak_m / (abs(total_m) + 1e-300);
    bool beyond = L > kp * r_bound;
    bool bad = !(RpL == RpL) || !(RmL == RmL);
    bool deep = depth_p > 1e12 || depth_m > 1e12;
    bool spur = beyond && (depth_p > 1e6 || depth_m > 1e6);
    if (bad || deep || spur) {
      unreliable[L] = 1;
      Rp[L] = 0.0;
      Rm[L] = 0.0;
    } else {
      unreliable[L] = 0;
      Rp[L] = RpL;
      Rm[L] = RmL;
    }
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerically stable evaluation of the retarded radial integrals for small photon wavenumber kp.
 *
 * For kp << 1 the closed-form branch decomposition of
 * dirac_rr_bf_retarded cancels terms of order kp^{-(m+1)} to leave
 * j_L(kp r) ~ (kp r)^L, a cancellation depth far beyond double precision
 * (for Z=1, E=1 eV, kp ~ 0.004 the depth reaches ~1e229, and Gamma(nu) hits
 * a pole as nu = gamma + s + i - m -> 0 for m ~ s + i).  Here j_L is instead
 * expanded in its convergent Taylor series,
 *  - j_L(x) = sum_t (-1)^t x^{L+2t} / (2^t t! (2L+2t+1)!!),
 * and the power-law integrals
 *  - Ip[p] = int_0^inf r^p P_c Q_b dr,
 *  - Im[p] = int_0^inf r^p Q_c P_b dr
 * are evaluated in closed form via the same Laplace-transform trick with the
 * positive parameter nu = gamma + s + i + p + 1 (never a pole, no
 * kp^{-(m+1)} factors):
 * - Ip[p] = pref_p sum_i b_coef[i] Re[Gamma(nu) mu^{-nu} ((1+c_f0) F1 + c_g0 F2)],
 * - Im[p] = pref_m sum_i a_coef[i] Re[Gamma(nu) mu^{-nu} (cQ_f F1 + cQ_g F2)],
 * with mu = lam - i kc, z2 = -2 i kc / mu, F1 = 2F1(a,nu;2s;z2), F2 = 2F1(a+1,nu;2s+1;z2).
 * The multipole integrals are then
 *  - Rp[L] = sum_t (-1)^t kp^{L+2t} / (2^t t! (2L+2t+1)!!) * Ip[L+2t],
 *  - Rm[L] = sum_t (-1)^t kp^{L+2t} / (2^t t! (2L+2t+1)!!) * Im[L+2t].
 *
 * Returns false if the t-series does not converge within tmax terms or
 * suffers deep cancellation -- the caller then keeps the closed-form result
 * instead.  Ip[p], Im[p] are cached per p, which is shared between the
 * (L, t) pairs with L + 2t = p.
 */
static bool dirac_rr_bf_retarded_taylor(double E_el_eV, double z, int kappa_c,
                                        const BoundWaveParams &bw, double kp,
                                        int Lmax, vector<double> &Rp,
                                        vector<double> &Rm) {
  ContinuumParams cp = dirac_continuum_params(E_el_eV, z, kappa_c);
  if (cp.k < 1e-30) {
    Rp.assign(Lmax + 1, 0.0);
    Rm.assign(Lmax + 1, 0.0);
    return true;
  }
  const int tmax = 40;
  if (bw.nr + Lmax + 2 * tmax > 1000) return false;  // coefficient arrays

  complex<double> ik(0.0, cp.k);
  complex<double> mu(bw.lam, -cp.k);  // lam - i kc
  complex<double> logmu = log(mu);
  complex<double> z2 = (-2.0 * ik) / mu;
  complex<double> c_f0 = (ik - cp.b1) / cp.c12;
  complex<double> c_g0 = -(ik * cp.a / cp.s) / cp.c12;
  complex<double> cQ_f = cp.q0 + cp.q2 * c_f0;
  complex<double> cQ_g = cp.q2 * c_g0;

  // log-factorials and log double factorials (2j+1)!! = prod_{k=0}^j (2k+1)
  vector<double> lfact(tmax + 1, 0.0);
  for (int k = 1; k <= tmax; k++) lfact[k] = lfact[k - 1] + log((double)k);
  vector<double> ldfact(Lmax + tmax + 1, 0.0);
  for (int j = 1; j <= Lmax + tmax; j++)
    ldfact[j] = ldfact[j - 1] + log(2.0 * j + 1.0);

  double logkp = log(kp);
  double s1p = cp.N * bw.C * exp(bw.gamma * log(bw.two_lam)) * sqrt(1.0 - bw.eps);
  double s1m = cp.N * bw.C * exp(bw.gamma * log(bw.two_lam)) * sqrt(1.0 + bw.eps);

  // Power integrals Ip[p], Im[p] (including the b_coef/a_coef sums and the
  // prefactors), cached lazily per p.  nu = gamma + s + i + p + 1 is always
  // positive, so no Gamma pole; the p actually reached stay well below the
  // Gamma-overflow threshold for the convergent t-sequences.
  int pmax = Lmax + 2 * tmax;
  vector<double> Ip(pmax + 1, 0.0), Im(pmax + 1, 0.0);
  vector<char> have(pmax + 1, 0);
  auto power_int = [&](int p) -> bool {
    if (have[p]) return true;
    complex<double> sp(0.0, 0.0), sm(0.0, 0.0);
    for (int i = 0; i <= bw.nr; i++) {
      double nu = bw.gamma + cp.s + i + p + 1.0;
      if (nu > 175.0) return false;  // Gamma(nu) would overflow
      double lgnu = std::lgamma(nu);
      complex<double> gp = exp(lgnu - nu * logmu);
      complex<double> F1, F2;
      if (!hypgeo2_cmplx(cp.a, complex<double>(nu, 0.0),
                         complex<double>(2.0 * cp.s, 0.0), z2, F1))
        return false;
      if (!hypgeo2_cmplx(cp.a + 1.0, complex<double>(nu, 0.0),
                         complex<double>(2.0 * cp.s + 1.0, 0.0), z2, F2))
        return false;
      sp += bw.b_coef[i] * gp * ((1.0 + c_f0) * F1 + c_g0 * F2);
      sm += bw.a_coef[i] * gp * (cQ_f * F1 + cQ_g * F2);
    }
    Ip[p] = s1p * sp.real();
    Im[p] = s1m * sm.real();
    have[p] = 1;
    return true;
  };

  Rp.assign(Lmax + 1, 0.0);
  Rm.assign(Lmax + 1, 0.0);
  for (int L = 0; L <= Lmax; L++) {
    double rp = 0.0, rm = 0.0;
    double ckp = 0.0, cmk = 0.0;  // Kahan compensations
    double absp = 0.0, absm = 0.0;
    int streak_p = 0, streak_m = 0;
    bool done = false;
    for (int t = 0; t <= tmax && !done; t++) {
      int p = L + 2 * t;
      // log |a_t|, a_t = (-1)^t kp^p / (2^t t! (2L+2t+1)!!)
      double loga = (double)p * logkp - (double)t * log(2.0) - lfact[t] -
                    ldfact[L + t];
      if (loga > 690.0) return false;
      double a = (t % 2) ? -exp(loga) : exp(loga);
      if (!power_int(p)) return false;
      double termp = a * Ip[p], termm = a * Im[p];
      absp += fabs(termp);
      absm += fabs(termm);
      {
        double y = termp - ckp;
        double tt = rp + y;
        ckp = (tt - rp) - y;
        rp = tt;
      }
      {
        double y = termm - cmk;
        double tt = rm + y;
        cmk = (tt - rm) - y;
        rm = tt;
      }
      streak_p = (fabs(termp) < 1e-14 * absp) ? streak_p + 1 : 0;
      streak_m = (fabs(termm) < 1e-14 * absm) ? streak_m + 1 : 0;
      if (streak_p >= 2 && streak_m >= 2) done = true;
    }
    if (!done) return false;  // t-series did not converge
    if (absp > 1e12 * fabs(rp) + 1e-300) return false;
    if (absm > 1e12 * fabs(rm) + 1e-300) return false;
    Rp[L] = rp;
    Rm[L] = rm;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical fallback for dirac_rr_bf_retarded
 *
 * The same Rp[L], Rm[L] integrals on a grid, with j_L(kp r) evaluated
 * pointwise.  Used whenever the closed-form m-branch sum flags a multipole
 * (NaN, cancellation depth beyond ~1e12, or a spurious moderate-cancellation
 * value past the turning point), and for every channel in builds without
 * FLINT.
 *
 * The continuum comes from the exact 1F1 construction
 * dirac_continuum_PQ_closed when HYDROCAL_HAVE_ACB_HYPGEOM is defined, and
 * from the ODE solver dirac_continuum_PQ otherwise.  The two agree in
 * r-dependence exactly -- they are proportional at every radius -- so the
 * only difference is the energy normalisation, and that is where the ODE
 * is weak: its envelope-matched normalisation stops converging at large
 * Sommerfeld parameter eta = Z(1+E/c^2)/k (high Z, low energy), leaving
 * (P_c,Q_c) uniformly too large by +0.72% at Z=92, E=1.5 eV (eta=277),
 * +0.44% at 10 eV, +0.24% at 100 eV, +0.12% at 1000 eV.  Since the cross
 * section is quadratic in the matrix element that becomes up to ~+1.4% in
 * sigma_RR, which is exactly the systematic excess over the Ichihara &
 * Eichler tables previously seen for Z=92 below 100 eV.  Builds without
 * FLINT still carry that bias (they have no exact 1F1 engine); it is
 * documented in sec:ichiharavalidation of DiracTransitionRates.tex.
 */
static void dirac_rr_bf_retarded_numeric(double E_el_eV, double z,
                                         int kappa_c,
                                         const BoundWaveParams &bw, double kp,
                                         int Lmax, vector<double> &Rp,
                                         vector<double> &Rm) {
  using hydroconst::alpha;
  using hydroconst::pi;
  double c = 1.0 / alpha;
  double E_el = E_el_eV / (2.0 * hydroconst::Ryd_eV);
  double k_grid = sqrt(E_el * (E_el + 2.0 * c * c)) / c;
  double r_max = max(50.0 / z, 60.0 / bw.lam);
  double k_max = max(k_grid, kp);
  int nper = (int)(k_max * r_max / (2.0 * pi)) + 1;
  int N = max(4000, 80 * nper + 1);
  double h = r_max / (N - 1);

  Rp.assign(Lmax + 1, 0.0);
  Rm.assign(Lmax + 1, 0.0);
  vector<double> jL(Lmax + 1);
  for (int ig = 0; ig < N; ig++) {
    double r = h * ig;
    double Pc = 0.0, Qc = 0.0, Pb = 0.0, Qb = 0.0;
    if (ig > 0) {
#ifdef HYDROCAL_HAVE_ACB_HYPGEOM
      // Exact 1F1 continuum: the ODE envelope-matched normalisation of
      // dirac_continuum_PQ is under-converged at large Sommerfeld parameter
      // (high Z, low energy) and biases every radial integral by the same
      // normalisation error, which enters the cross section twice.
      dirac_continuum_PQ_closed(E_el_eV, z, kappa_c, r, Pc, Qc);
#else
      dirac_continuum_PQ(E_el_eV, z, kappa_c, r, Pc, Qc);
#endif
      double rpow = pow(bw.two_lam * r, bw.gamma);
      double expfac = exp(-bw.lam * r);
      double Fsum = 0.0, Gsum = 0.0;
      for (int i = 0; i <= bw.nr; i++) {
        Fsum += bw.a_coef[i] * pow(r, i);
        Gsum += bw.b_coef[i] * pow(r, i);
      }
      Pb = bw.C * sqrt(1.0 + bw.eps) * rpow * expfac * Fsum;
      Qb = bw.C * sqrt(1.0 - bw.eps) * rpow * expfac * Gsum;
    }
    for (int L = 0; L <= Lmax; L++) jL[L] = sbesj(L, kp * r);
    double w = (ig == 0 || ig == N - 1) ? 0.5 * h : h;
    for (int L = 0; L <= Lmax; L++) {
      Rp[L] += w * Pc * Qb * jL[L];
      Rm[L] += w * Qc * Pb * jL[L];
    }
  }
}

////////////////////////////////////////////////////////////////////////////
/**
 * Per-channel fully-retarded matrix element sum
 *
 * tot(kappa_c) = sum_{ms2a, ms2b, pol} |M|^2
 * assembled exactly as the bound-bound retarded_total in DiracRate.cxx: the
 * phase conventions for Mx, My are
 * - Mx = i[(A0B3) + (A1B2) - (A2B1) - (A3B0)]
 * - My =    [(A0B3) - (A1B2) - (A2B1) + (A3B0)]
 * with (large a x small b: +i) and (small a x large b: -i).
 */
 static double retarded_channel_tot(double E_el_eV, double z, int kappa,
                                   int kappa_c, const BoundWaveParams &bw,
                                   double kp, bool numeric, bool &ok) {
  int lc, lb;
  double jc, jb;
  lj_from_kappa(kappa_c, lc, jc);
  lj_from_kappa(kappa, lb, jb);
  auto l_of = [](int kp_in) {
    int l;
    double j;
    lj_from_kappa(kp_in, l, j);
    return l;
  };
  int Lmax = max(l_of(kappa_c) + l_of(-kappa), l_of(-kappa_c) + lb);

  vector<double> Rp, Rm;
  if (numeric) {
    dirac_rr_bf_retarded_numeric(E_el_eV, z, kappa_c, bw, kp, Lmax, Rp, Rm);
  } else {
    // The branch form is well-conditioned for kp of order unity and above but
    // cancels to garbage below (kp < ~0.5); the Taylor form is well-condi-
    // tioned exactly there.  Below the threshold evaluate the whole channel
    // in the Taylor series; otherwise take the closed form and patch only
    // the multipoles the branch form flags as unreliable.
    if (kp < 0.5) {
      if (!dirac_rr_bf_retarded_taylor(E_el_eV, z, kappa_c, bw, kp, Lmax, Rp, Rm)) {
        ok = false;
        return 0.0;
      }
    } else {
      vector<char> unreliable;
      if (!dirac_rr_bf_retarded(E_el_eV, z, kappa_c, bw, kp, Lmax, Rp, Rm, unreliable)) {
        ok = false;
        return 0.0;
      }
      // Multipoles the branch form flagged (NaN, cancellation depth beyond
      // ~1e12, or spurious moderate-cancellation values beyond the turning
      // point) are recomputed on the grid.  The direct quadrature has no
      // branch-conditioning at all and reproduces the validated channel
      // totals to well below 1%; it costs a small fraction of a second per
      // channel (the grid is shared through the pointwise j_L loop), so the
      // fallback is used freely whenever the closed form is not trustworthy.
      bool need_numeric = false;
      for (size_t L = 0; L < unreliable.size(); L++)
        if (unreliable[L]) need_numeric = true;
      if (need_numeric) {
        dirac_rr_bf_retarded_numeric(E_el_eV, z, kappa_c, bw, kp, Lmax, Rp, Rm);
      }
    }
  }

  double tot = 0.0;
  int ma_min = -(int)(2.0 * jc), ma_max = (int)(2.0 * jc);
  int mb_min = -(int)(2.0 * jb), mb_max = (int)(2.0 * jb);
  for (int ms2a = ma_min; ms2a <= ma_max; ms2a += 2) {
    AngComp A0, A1, A2, A3;
    spinor_harm(kappa_c, ms2a, A0, A1);
    spinor_harm(-kappa_c, ms2a, A2, A3);
    for (int ms2b = mb_min; ms2b <= mb_max; ms2b += 2) {
      AngComp B0, B1, B2, B3;
      spinor_harm(kappa, ms2b, B0, B1);
      spinor_harm(-kappa, ms2b, B2, B3);
      complex<double> S0(0.0, 0.0), S1(0.0, 0.0), S2(0.0, 0.0), S3(0.0, 0.0);
      {
        vector<complex<double> > ac = ang_int(A0.l, A0.m, B3.l, B3.m);
        for (size_t L = 0; L < ac.size(); L++)
          if (ac[L] != complex<double>(0.0, 0.0)) S0 += ac[L] * Rp[L];
      }
      {
        vector<complex<double> > ac = ang_int(A1.l, A1.m, B2.l, B2.m);
        for (size_t L = 0; L < ac.size(); L++)
          if (ac[L] != complex<double>(0.0, 0.0)) S1 += ac[L] * Rp[L];
      }
      {
        vector<complex<double> > ac = ang_int(A2.l, A2.m, B1.l, B1.m);
        for (size_t L = 0; L < ac.size(); L++)
          if (ac[L] != complex<double>(0.0, 0.0)) S2 += ac[L] * Rm[L];
      }
      {
        vector<complex<double> > ac = ang_int(A3.l, A3.m, B0.l, B0.m);
        for (size_t L = 0; L < ac.size(); L++)
          if (ac[L] != complex<double>(0.0, 0.0)) S3 += ac[L] * Rm[L];
      }
      complex<double> c03 = A0.coeff * B3.coeff * S0;
      complex<double> c12 = A1.coeff * B2.coeff * S1;
      complex<double> c21 = A2.coeff * B1.coeff * S2;
      complex<double> c30 = A3.coeff * B0.coeff * S3;
      complex<double> Mx = complex<double>(0.0, 1.0) * (c03 + c12 - c21 - c30);
      complex<double> My = c03 - c12 - c21 + c30;
      tot += norm(Mx) + norm(My);
    }
  }
  ok = true;
  return tot;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief relativistic hydrogenic RR cross section (fully retarded)
 *
 * Fully-retarded (plane-wave, all multipoles) RR cross section per vacancy
 * in cm^2. The retarded operator e^{-ik_ph r cosT} couples all multipoles
 * and (unlike the E1 dipole cross seection) both parities, and is the
 * physically complete description; it matches the "exact relativistic" tables
 * of Ichihara & Eichler (which include all multipole orders and finite nuclear
 * size) to the level of the point-nucleus idealisation.  
 *
 * The closed-form retarded radial integrals require the exact 2F1 engine from
 * FLINT; without it (HYDROCAL_HAVE_ACB_HYPGEOM undefined) every channel is
 * evaluated with the grid quadrature dirac_rr_bf_retarded_numeric, which
 * needs no hypergeometric engine at all (slower, but validated against the
 * closed form to well below 1%).  Without FLINT that quadrature has to take
 * its continuum from the ODE solver, whose energy normalisation is
 * under-converged at large Sommerfeld parameter: at Z=92 it overestimates
 * sigma_RR by ~1.4% at 1.5 eV, ~0.5% at 100 eV and ~0.2% at 1 keV (the bias
 * vanishes as the energy rises).  With FLINT the quadrature uses the exact
 * 1F1 continuum and agrees with the Ichihara & Eichler tables to within their
 * own three-digit rounding over 1.5 eV -- 10 MeV.
 *
 * The matrix element sum runs over the complete continuum channel
 * set |j_c - j_b| = 0, 1, 2, ... (each j_c realised by both
 * kappa_c = +-(j_c + 1/2)); the per-channel contributions are
 * positive-definite, so the shells are added until one shell adds less than
 * rel_tol of the running total.  Cross section via
 *   sigma_RR = [2 pi^2/(alpha omega)] * sum_tot/(2j_b+1)
 *              * [alpha^2 omega^2/(2 k_el^2)]   [a0^2]
 * (the fully-retarded analogue of the long-wavelength Milne formula of
 * dirac_rr_xsec_dipole; the alpha omega normalisation reproduces the dipole result
 * as k_ph -> 0).
 *
 * @param E_el_eV electron kinetic energy (eV)
 * @param z       nuclear charge
 * @param n       principal quantum number of the captured state
 * @param kappa   relativistic angular quantum number of the captured state
 *
 * @return sigma_RR in cm^2, never negative
 * @return  0.0 for a non-existent state, non-positive energy, or an
 *          evaluation that could not be completed
 */
double dirac_rr_xsec_retarded(double E_el_eV, double z, int n, int kappa)
{
  using hydroconst::alpha;
  using hydroconst::a0_cm;
  using hydroconst::Ryd_eV;
  using hydroconst::pi;

  if (E_el_eV <= 0.0) return 0.0;
  // Reject non-existent bound states (l >= n)
  if (!valid_kappa(n, kappa)) return 0.0;
  int lb; double jb; // bound state angular momenta
  lj_from_kappa(kappa, lb, jb); 

  double E_el = E_el_eV / (2.0 * Ryd_eV); // Hartree
  double E_bind = dirac_binding_energy(z, n, kappa); // eV, negative
  double omega_eV = E_el_eV + fabs(E_bind);
  double omega = omega_eV / (2.0 * Ryd_eV); // Hartree
  double kp = omega * alpha;  // photon wavenumber in a.u. (omega Hartree)
  double k_el = sqrt(E_el * (E_el + 2.0 / (alpha * alpha))) * alpha;

  // Bound-state wavefunction parameters
  BoundWaveParams bw = dirac_bound_wave_params(z, n, kappa);
  double za = z * alpha;

  // Validation cross-check: HYDROCAL_DIRAC_RETARDED_NUMERIC=1 forces the
  // grid quadrature (dirac_rr_bf_retarded_numeric) for every channel.
  // Without FLINT the closed forms (dirac_rr_bf_retarded, _taylor) cannot be
  // trusted -- they sit on hypgeo2_cmplx, whose portable engine may return
  // false or an untrustworthy value in exactly the deep-cancellation regime
  // the radial integrals live in -- so the quadrature, which needs no
  // hypergeometric engine, becomes the default instead of an error return.
  static const bool numeric = [] {
    const char *env = getenv("HYDROCAL_DIRAC_RETARDED_NUMERIC");
    if (env != nullptr && *env != '\0' && strcmp(env, "0") != 0) return true;
#ifndef HYDROCAL_HAVE_ACB_HYPGEOM
    return true;  // no exact 2F1 engine: grid quadrature for every channel
#else
    return false;
#endif
  }();
  static const bool dbg = [] {
    const char *env = getenv("HYDROCAL_DIRAC_DEBUG");
    return env != nullptr && *env != '\0' && strcmp(env, "0") != 0;
  }();

  const double rel_tol = 1e-4;
  const int dj_max = 99;
  // Per-shell channel sums are independent, so the expensive part (the
  // closed-form 2F1 evaluations, one per channel) is parallelised over dj.
  // The shells are then accumulated in ascending order, preserving the
  // sequential convergence semantics.
  vector<double> shell_arr(dj_max + 1, 0.0);
  vector<int> any_arr(dj_max + 1, 0);
#pragma omp parallel for schedule(dynamic, 4)
  for (int dj = 0; dj <= dj_max; dj++) {
    double shell = 0.0;
    int any = 0;
    for (int sg = -1; sg <= 1; sg += 2) {
      if (dj == 0 && sg < 0) continue;  // jc = jb counted once at dj = 0
      double jc = jb + sg * dj;
      if (jc < 0.5) continue;
      int jh = (int)(2.0 * (jc + 0.5) + 0.5);  // 2*jc+1 (exact integer)
      int kapc = jh / 2;                        // jc + 1/2
      for (int ks = -1; ks <= 1; ks += 2) {
        int kappa_c = ks * kapc;
        if (kappa_c == 0) continue;
        double gamma_c = sqrt(kappa_c * kappa_c - za * za);
        if (gamma_c <= 0.0) continue;
        bool ok = false;
        double ch = retarded_channel_tot(E_el_eV, z, kappa, kappa_c, bw, kp, numeric, ok);
        if (dbg)
          fprintf(stderr, "  dj=%d jc=%.2f kc=%d ch=%.6e ok=%d\n", dj, jc,
                  kappa_c, ch, (int)ok);
        if (ok) {
          any = 1;
          shell += ch;
        }
      }
    }
    shell_arr[dj] = shell;
    any_arr[dj] = any;
    if (dbg) fprintf(stderr, "dj=%d shell=%.6e\n", dj, shell);
  }
  double tot = 0.0;
  bool any_ok = false;
  bool converged = false;
  for (int dj = 0; dj <= dj_max; dj++) {
    if (any_arr[dj]) any_ok = true;
    tot += shell_arr[dj];
    if (any_ok && dj > 0 && shell_arr[dj] < rel_tol * tot) {
      converged = true;
      break;
    }
  }
  if (!any_ok || tot < 1e-300) return 0.0;
  if (!converged)
    fprintf(stderr,
            "WARNING: retarded RR channel series not converged at E=%.4g eV "
            "Z=%g n=%d kappa=%d (dj cap %d reached)\n",
            E_el_eV, z, n, kappa, dj_max);

  double sigma_au = 2.0 * pi * pi / (alpha * omega) * tot / (2.0 * jb + 1.0);
  sigma_au *= alpha * alpha * omega * omega / (2.0 * k_el * k_el);
  return sigma_au * a0_cm * a0_cm;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief relativistic hydrogenic RR cross section (dipole approximation)
 *
 * This is a long-wavelength (dipole, E1-only) cross section. 
 * Uses the Milne (detailed-balance relation:
 *   sigma_RR = pi^2 alpha^3 a0^2 x (2j+1) x (omega/k)^2 x (df/dE)
 * where df/dE is the bound-free oscillator strength per unit energy
 * computed from the Dirac electric-dipole matrix elements, and k is the
 * relativistic electron momentum.
 *
 * The radial element R_bf = int r (P_b P_c + Q_b Q_c) dr omits the retarded
 * multipole operators (j_L(omega r/c) factors) and the higher multipoles
 * (E2, M1, ...), so the result is reliable only where omega r/c << 1 over
 * the capture volume, i.e. for photon energies well below ~c/z (in
 * atomic units) times the inverse size of the bound state.  Above that
 * the tabulated exact relativistic cross sections of Ichihara & Eichler
 * (see sigma_IchiharaEichlerRR) grow relative to this dipole value as
 * retardation and the higher multipoles set in; for example at Z=110,
 * n=3, 2 MeV this dipole cross section is ~2e-3 of the exact value.
 *
 *
 * @param E_el_eV electron kinetic energy (eV)
 * @param z       nuclear charge
 * @param n       principal quantum number of the captured state
 * @param kappa   relativistic angular quantum number of the captured state
 * @param retarded_flag fully retarded (if true), dipole approximation (if false)
 *
 * @return sigma_RR in cm^2, 0 for a non-existent state or non-positive
 *         energy
 */
double dirac_rr_xsec_dipole(double E_el_eV, double z, int n, int kappa) {
  using hydroconst::alpha;
  using hydroconst::a0_cm;
  using hydroconst::Ryd_eV;

  if (E_el_eV <= 0.0) return 0.0;
  // Reject non-existent bound states (l >= n)
  if (!valid_kappa(n, kappa)) return 0.0;
  double E_el = E_el_eV / (2.0 * Ryd_eV); // Hartree

  // Continuum-wavefunction selection (hybrid by default).
  //
  // Two continuum solvers are available: the numerical ODE integration of
  // dirac_continuum_PQ (with its envelope-matched energy normalisation) and
  // the exact closed-form 1F1 construction of dirac_continuum_PQ_closed
  // (sec:closedform of DiracTransitionRates.tex).  Neither is best
  // everywhere:
  //   * The ODE energy normalisation is under-converged at large Sommerfeld
  //     parameter eta = z(1+E_el/c^2)/k (low E_kin, high z): its
  //     envelope-convergence loop stops while the RMS amplitude is still
  //     drifting, leaving up to ~40% normalisation error at Z=110,
  //     E_kin=1e-5 eV (sec:ichiharavalidation of DiracTransitionRates.tex).
  //   * The closed form is exact, but its double/quad-precision 1F1
  //     evaluation (~1e-8 relative) cannot resolve the deep-cancellation
  //     bound-free integrals of the high-energy regime, where the ODE plus
  //     Gordon fallback (sec:rrxsec) is retained.
  // Default (env unset): hybrid -- closed form for eta > closed_eta_thr
  // (the low-energy/high-z corner, where it is exact), ODE otherwise.
  // HYDROCAL_DIRAC_CLOSED=0 forces the ODE everywhere (legacy default);
  // any other non-empty value forces the closed form everywhere (validation
  // cross-check; slower per grid point, and the exact 1F1 engine is used
  // only when built with FLINT).
  static const double closed_eta_thr = 350.0;
  static const int closed_mode = [] {
    const char *env = getenv("HYDROCAL_DIRAC_CLOSED");
    if (env == nullptr || *env == '\0') return -1;  // hybrid
    if (strcmp(env, "0") == 0) return 0;            // force ODE
    return 1;                                       // force closed
  }();
  double k_sel = sqrt(E_el * (E_el + 2.0 / (alpha * alpha))) * alpha;
  double eta = z * (1.0 + E_el * alpha * alpha) / k_sel;
  bool use_closed_form = (closed_mode == 1) ||
                         (closed_mode == -1 && eta > closed_eta_thr);

  int l;
  double j;
  lj_from_kappa(kappa, l, j);

  double E_bind = dirac_binding_energy(z, n, kappa); // eV, negative
  double omega_eV = E_el_eV + fabs(E_bind);
  double omega = omega_eV / (2.0 * Ryd_eV); // Hartree

  // Bound-state wavefunction parameters
  BoundWaveParams bw = dirac_bound_wave_params(z, n, kappa);
  double za = z * alpha;
  //
  // Two ways of computing the bound-free element
  //   R_bf = int_0^inf r [P_b P_c + Q_b Q_c] dr
  // are available.  The analytic (relativistic-Gordon) term-sum of
  // dirac_rr_bf_analytic evaluates R_bf in closed form through the
  // Laplace transform of 1F1 (a finite sum of 2F1 values), with no
  // quadrature and no deep-cancellation loss at any energy; it needs the
  // exact 2F1 engine (FLINT) and is the default when that is available.
  // The legacy path integrates the wavefunctions on a logarithmic grid
  // (trapezoidal rule), substituting the nonrelativistic Gordon element
  // with a runtime-calibrated conversion factor where the quadrature
  // cancels to below its roundoff; it is retained as the fallback for
  // builds without FLINT and as a cross-check.
  // HYDROCAL_DIRAC_GORDON=0 forces the legacy path; =1 forces the analytic
  // path even without FLINT (portable 2F1, best effort).
  static const int gordon_mode = [] {
    const char *env = getenv("HYDROCAL_DIRAC_GORDON");
    if (env == nullptr || *env == '\0') {
#ifdef HYDROCAL_HAVE_ACB_HYPGEOM
      return 1;  // analytic by default when FLINT is available
#else
      return 0;  // legacy by default otherwise
#endif
    }
    return (strcmp(env, "0") == 0) ? 0 : 1;
  }();

  // Sum over E1-coupled continuum channels: kappa_c = kappa +/- 1 and
  // kappa_c = -kappa (e.g. s1/2 <-> p1/2)
  double total_me_sq = 0.0;
  int kappa_candidates[] = {kappa + 1, kappa - 1, -kappa};

  for (int ik = 0; ik < 3; ik++) {
    int kappa_c = kappa_candidates[ik];
    // Check validity
    int lc; double jc;
    lj_from_kappa(kappa_c, lc, jc);
    // Parity selection: l + lc + 1 must be even
    if ((l + lc + 1) % 2 != 0) continue;
    // Finite angular momentum: |jc - j| <= 1, j > 0
    if (fabs(jc - j) > 1.01) continue;
    if (j < 0.01 && jc < 0.01) continue;

    if (kappa_c == 0) continue;
    double gamma_c = sqrt(kappa_c * kappa_c - za * za);
    if (gamma_c <= 0.0) continue;

    // 3j symbol for the angular factor
    double tj = ThreeJ(jc, 1.0, j, -0.5, 0.0, 0.5);
    if (fabs(tj) < 1e-15) continue;

    // Bound-free radial integral.  Analytic (relativistic-Gordon) term-sum
    // by default (see the gordon_mode selection above): a finite sum of
    // 2F1 values, exact at every energy and free of the deep-cancellation
    // loss of the numerical quadrature.  The legacy numerical path is the
    // fallback when the analytic engine is unavailable or did not produce
    // a trustworthy 2F1 (portable engine only).
    double R_bf = 0.0;
    if (gordon_mode == 1 &&
        dirac_rr_bf_analytic(E_el_eV, z, kappa_c, bw, R_bf)) {
      // Analytic path: nothing more to do.
    } else { // +++++++++++++++++++ start of numerical path  ++++++++++++++++++++++++++++++++ 
      // Legacy numerical path: integrate the wavefunctions on a
      // logarithmic grid (trapezoidal rule), substituting the
      // nonrelativistic Gordon element with a runtime-calibrated
      // conversion factor where the quadrature cancels below its roundoff.

      // Outer radial cutoff for the bound-free overlap integral: the bound
      // state decays as exp(-lam*r) (times a degree-nr polynomial), so its
      // extent grows with n (lam ~ z/n). A cutoff fixed at 50/z regardless
      // of n truncates that tail before it has decayed away for n>=2 --
      // leaving the integral incomplete right where the true R_bf is a
      // near-total cancellation between the bound state's oscillatory
      // overlap with the continuum, and an incomplete cancellation reads
      // as a large spurious residual. Require enough e-foldings (60) that
      // the tail is negligible to well below double precision at r_max,
      // whatever n is.
      double r_max = std::max(50.0 / z, 60.0 / bw.lam);

      // Exact (Gordon) bound-free radial matrix elements as an analytic
      // fallback. For hydrogenic ions the bound-free integral can be
      // evaluated in closed form; here we compute the nonrelativistic
      // Gordon element via CalcRbc at the scaled energy e = E/(z^2 Ry) and
      // rescale it to this code's R_bf convention with
      // R_bf = c*(2/z^2)*R_gordon, where the residual factor c is
      // calibrated at runtime from channels whose numerical integration is
      // reliable (see below). The relativistic correction to the radial
      // element is O((Z*alpha)^2) ~ 5e-5 for Z=1, so the substitution is
      // asymptotically correct in the high-energy regime where the
      // trapezoidal R_bf suffers catastrophic cancellation.
      double e_scaled = E_el_eV / (z * z * Ryd_eV);
      vector<double> Rbcp(n + 1, 0.0), Rbcm(n + 1, 0.0);
      if (e_scaled > 0.0) CalcRbc(1.0 / sqrt(e_scaled), n, Rbcp, Rbcm);

      // Runtime calibration of the conversion factor
      // c = R_bf/((2/z^2)*R_gordon), averaged over all reliable channels
      // computed at the current (E, z). The scatter of c across channels
      // measures whether a single scalar can represent the relativistic
      // correction at all (it can for small (Z*alpha)^2, but not for large
      // Z where the correction is strongly channel-dependent); see the
      // substitution guard below.
      static double cal_E = -1.0, cal_z = -1.0, cal_c_sum = 0.0;
      static double cal_c_sqsum = 0.0;
      static int cal_c_n = 0;
      if (cal_E != E_el_eV || cal_z != z) {
        cal_E = E_el_eV;
        cal_z = z;
        cal_c_sum = 0.0;
        cal_c_sqsum = 0.0;
        cal_c_n = 0;
      }

      // The Gordon substitution is only trustworthy where the
      // non-relativistic element approximates the Dirac element: small
      // (Z*alpha)^2, and a calibrated c that is consistent across reliable
      // channels (relative scatter < 20%). Outside this regime an
      // unreliable channel is dropped rather than silently assigned a wrong
      // value. (Z*alpha)^2 > 0.1 is Z > ~43; the fully relativistic element
      // has no closed form here.
      double za2 = za * za;
      bool gordon_ok = za2 < 0.1;

      // Radial grid: logarithmic from r_min to r_max. The continuum
      // wavefunction oscillates with a *fixed absolute* period 2*pi/k_c,
      // but a log grid's absolute point spacing grows proportionally with r
      // (dr = r * ln(r_max/r_min)/(ng-1)) -- so a point count sized only as
      // 32*r_max/period_grid (as if spacing were uniform) under-resolves
      // the oscillation everywhere except near r_min by a factor of
      // ln(r_max/r_min) (>15 for the r_min/r_max used here). That aliases
      // the trapezoidal R_bf integral specifically where the bound state
      // still has weight out at large r -- negligible for compact n=1
      // states, but large and energy-dependent (hence the erratic ratios)
      // for n=3 states whose radial extent reaches well into the
      // coarsely-sampled outer grid. Include the log-Jacobian factor so the
      // worst-case (largest r) spacing is still resolved.
      double r_min = 1e-6;
      double c_grid = 1.0 / alpha;
      double k_grid = sqrt(E_el * (E_el + 2.0 * c_grid * c_grid)) / c_grid;
      double period_grid =
          k_grid > 0.0 ? 2.0 * hydroconst::pi / k_grid : r_max;
      double log_range = log(r_max / r_min);
      int ng = std::max(2000, (int)(32.0 * r_max * log_range / period_grid) + 1);
      ng = std::min(ng, 2000000);

      // Precompute bound state on grid
      vector<double> P_b(ng), Q_b(ng), P_c(ng), Q_c(ng), rr(ng);
      for (int ig = 0; ig < ng; ig++) { 
        double t = (double)ig / (ng - 1);
        double r_grid = r_min * exp(t * log(r_max / r_min));
        rr[ig] = r_grid;

        // Bound state
        double rpow = pow(bw.two_lam * r_grid, bw.gamma);
        double expfac = exp(-bw.lam * r_grid);
        double Fsum = 0.0, Gsum = 0.0;
        for (int i = 0; i <= bw.nr; i++) {
          Fsum += bw.a_coef[i] * pow(r_grid, i);
          Gsum += bw.b_coef[i] * pow(r_grid, i);
        }
        P_b[ig] = bw.C * sqrt(1.0 + bw.eps) * rpow * expfac * Fsum;
        Q_b[ig] = bw.C * sqrt(1.0 - bw.eps) * rpow * expfac * Gsum;

        // Continuum state
        if (use_closed_form)
          dirac_continuum_PQ_closed(E_el_eV, z, kappa_c, r_grid, P_c[ig],
                                    Q_c[ig]);
        else
          dirac_continuum_PQ(E_el_eV, z, kappa_c, r_grid, P_c[ig], Q_c[ig]);
      }  // end (for ig ...)

      // Integrate ∫ r (P_b P_c + Q_b Q_c) dr using trapezoidal rule
      R_bf = 0.0;
      double S_max = 0.0;  // peak |running partial sum| (cancellation probe)
      for (int ig = 1; ig < ng; ig++) {
        double dr = rr[ig] - rr[ig-1];
        double f_prev = rr[ig-1] * (P_b[ig-1] * P_c[ig-1] +
                                    Q_b[ig-1] * Q_c[ig-1]);
        double f_cur = rr[ig] * (P_b[ig] * P_c[ig] + Q_b[ig] * Q_c[ig]);
        R_bf += 0.5 * dr * (f_prev + f_cur);
        S_max = std::max(S_max, fabs(R_bf));
      } // end for(ig ...)

      // Reliability test: when the oscillatory bound-continuum overlap
      // cancels almost completely (|R_bf| << peak |partial sum|), the
      // trapezoidal result is dominated by accumulated roundoff and must
      // not be trusted. Substituting the exact analytic (Gordon) element
      // via CalcRbc at the scaled energy e = E/(z^2 Ry) is then
      // asymptotically correct to O((Z*alpha)^2). Mapping (verified
      // numerically): the l -> l+1 channel uses Rbcm[l+1], the l -> l-1
      // channel uses Rbcp[l].
      bool reliable = (S_max > 0.0) && (fabs(R_bf) > 1e-8 * S_max);
      double g_el = 0.0;
      if (lc == l + 1) g_el = Rbcm[l + 1];
      else if (lc == l - 1) g_el = Rbcp[l];
      double R_bf_gordon = (2.0 / (z * z)) * g_el;
      if (reliable) {
        // Calibrate the residual factor c from reliable channels of this
        // (E, z) (c drifts with the relativistic kinematics, ~1.0003 at
        // low energy to ~0.978 at 10 keV for Z=1)
        if (R_bf_gordon != 0.0) {
          double c = R_bf / R_bf_gordon;
          cal_c_sum += c;
          cal_c_sqsum += c * c;
          cal_c_n++;
        }
      } else if (R_bf_gordon != 0.0) {
        // Substitution guard: apply the Gordon fallback only where it is a
        // legitimate asymptotic approximation -- small (Z*alpha)^2 -- and
        // where the calibrated c is consistent across the reliable channels
        // (relative scatter < 20%). With no calibration samples the pure
        // Gordon element (c = 1) is still correct to O((Z*alpha)^2). When
        // the guard fails (e.g. Z=92 at ~1 MeV, where the per-channel
        // relativistic restructuring scatters c from -0.5 to +4.3), the
        // unreliable channel is dropped rather than silently assigned a
        // wrong value.
        double c_cal = 1.0;
        bool scatter_ok = true;
        if (cal_c_n > 0) {
          double c_mean = cal_c_sum / cal_c_n;
          c_cal = c_mean;
          double var = cal_c_sqsum / cal_c_n - c_mean * c_mean;
          if (var < 0.0) var = 0.0;
          scatter_ok = (c_mean > 0.0) && (sqrt(var) / c_mean < 0.2);
        }
        if (gordon_ok && scatter_ok) R_bf = c_cal * R_bf_gordon;
        else R_bf = 0.0;
      } // end if(R_bf_gordon!=)
    } // +++++++++++++++++++ end of numerical path  ++++++++++++++++++++++++++++++++++ 

    // Reduced matrix element squared:
    // |⟨κ_c||C^(1)||κ⟩|² = (2jc+1)(2j+1) × [3j(jc,1,j;-1/2,0,1/2)]²
    double ang = (2.0 * jc + 1.0) * (2.0 * j + 1.0) * tj * tj;

    total_me_sq += ang * R_bf * R_bf;
  }

  // total_me_sq is a sum of non-negative channel contributions ang*R_bf^2,
  // each of which is either a reliable numeric integral or the exact Gordon
  // fallback, so a very low floor is safe. Higher-l high-energy channels
  // legitimately reach R_bf^2 well below 1e-30 and must not be discarded.
  if (total_me_sq < 1e-100) return 0.0;

  // Bound-free oscillator strength per unit energy (Hartree⁻¹):
  //   df/dE = (2/3) × ω × |⟨κ_f||C^(1)||κ_i⟩|² × R² / (2j_i+1)
  double dfdE = (2.0 / 3.0) * omega * total_me_sq / (2.0 * j + 1.0);

  // Milne (detailed-balance) relation between photoionisation and
  // radiative recombination:
  //   σ_PI(ω)      = 2π²α a₀² × (df/dE)                        [a₀²]
  //   σ_RR(E_el)   = σ_PI(ω) × (g_bound / g_continuum) × (ω/(p c))²
  // with g_continuum = 2 (free-electron spin, ion core g_+ = 1), and p the
  // *relativistic* electron momentum, p c = ħk with k² = E_el(E_el + 2c²)/c²
  // (same k used by the continuum solver above). df/dE above is *already*
  // per one bound magnetic substate (the explicit (2j+1) in the reduced
  // matrix element ang = (2jc+1)(2j+1) tj² is divided back out just above),
  // so g_bound here is 1 for this per-vacancy cross section -- not (2j+1),
  // which would double count the degeneracy already folded into df/dE.
  // Combining:
  //   σ_RR = π² α³ a₀² × ω² / k² × (df/dE)
  // which reduces, for k² -> 2E_el (E_el << c²), to the standard
  // nonrelativistic Milne cross section per vacancy.
  double c2 = 1.0 / (alpha * alpha);
  double k2 = E_el * (E_el + 2.0 * c2) / c2;
  double sigma_au = hydroconst::pi * hydroconst::pi * alpha * alpha * alpha *
                    a0_cm * a0_cm * omega * omega / k2 * dfdE;

  // Convert from a.u. to cm²: 1 a₀² = a0_cm²
  // The formula above already uses a0_cm, so sigma_au is in cm²
  if (sigma_au < 0.0) return 0.0;
  return sigma_au;
}
