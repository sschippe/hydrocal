// $Id: sigmarr.h 602 2019-03-28 17:47:11Z iamp $
// Calculation of hydrogenic radiative recombination cross sections 
// semiclassically and quantum mechanically. 

// Copyright: 1997,2005 by Stefan.E.Schippers@iamp.physik.uni-giessen.de 



#include "fele.h"


/** this type definition allows functions returning a cross section to be passed as parameters
 */
typedef double (*SIGMARR)(double, double, double* , int, int, int, int);

/** convolution of a delta peak like cross section
 *
 * @return alpha_DR(E)
 */
double deltapeak(FELE fele, double er, double eV, double ktpar, double ktperp, double a);

//////////////////////////////////////////////////////////////////////////
/**
 * @brief initializes interpolation arrays
 */
void init_sigma_interpolation(double *energy, double *xsec, int npts);

//////////////////////////////////////////////////////////////////////////
/**
 * @brief cleans up interpolation arrays
 */
void delete_sigma_interpolation();

//////////////////////////////////////////////////////////////////////////
/**
 * @brief returns interpolated cross section
 */
double sigma_interpolated(double e, double z, double *fraction, 
			  int use_fraction, int nmin, int nmax, int lmin);


/** sigma_DR(E) ~ Lorentzian or ~ Lorentzian/E depending on the definition of the parameter fraction[3]. 
 *
 * @return E*sigma_DR(E)
 */
double sigmalorentzian(double e, double er, double *fraction,
                   int use_fraction, int nmin, int nmax, int lmin);

/** test cross section yielding alpharr = 1.0
 */
double sigmarrtest(double e, double z, double *fraction, 
                   int use_fraction, int nmin, int nmax, int lmin);


/** lim e->0 e*sigma (eV cm^2)
 */
double sigmarrqm(double z, int n, int l);

/**  e*sigma (eV cm^2)
 */
double sigmarrqm(double eV, double z, int n, int l);

/**  hydrogenic radiative
 *   recombination cross section times energy using bound free oscillator 
 *   strengths at energy eV, with an effective charge z. The factors fraction 
 *   account for a partial filling of the lowest free shell denoted by nmin as 
 *   well as for the cutoff of Rydberg states in the dipole magnets.
 *   nmax is the cutoff quantum number.
*/
double sigmarrqm(double eV, double z, double* fraction, 
                 int use_fraction, int nmin, int nmax, int lmin);

/**  is the semi-classical Bethe Salpeter
 *  Formula for the radiative recombination cross section times energy 
 *   with Stobbe corrections at energy eV, with an effective charge z. The 
 *   factors fraction account for a partial filling of the lowest free shell 
 *   denoted by nmin as well as for the cutoff of Rydberg states in the dipole 
 *   magnets.  nmax is the cutoff quantum number.
*/
double sigmarrscl(double eV, double z, double* fraction, int use_fraction, 
                  int nmin, int nmax, int lmin);

/**  sigmarrscl up to n2-1, from n2 to nmax approximation by integral over n
*/
double sigmarrscl2(double eV, double z, double *fraction, 
                   int n2, int nmin, int nmax, int lmin);



/**  is the quantum mechanical hydrogenic 
 *  radiative recombination cross section times energy at energy eV, with an 
 *   effective charge z. The factors fraction account for a partial filling of 
 *   the lowest free shell denoted by nmin as well as for the cutoff of Rydberg 
 *   states in the dipole magnets.  nmax is the cutoff quantum number.
 *   For n>45 the quantum mechanical cross section is calculated by multiplying
 *   the semiclassical cross section with the energy dependend Gaunt factor.
 *   Use is made of the fact that the Gaunt factor becomes nearly independent
 *   from n at large n by using the Gaunt factor for n=50 also for higher
 *   n values. This approximates the correct quantum mechanical result to better
 *   than 1% at all energies.
*/
double sigmarraqm(double eV, double z, double* fraction, int use_fraction,
                  int nmin, int nmax, int lmin);

/**  calculates the Stobbe correction factor to the semiclassical
 *  radiative recombination cross section into the hydrogenic level with
 *   main quantum number n. 
*/
double stobbe(int n);

void calc_stobbe(void);

/** 
 * @brief Semiclassical RR rate coeffcient in a plasma
 * 
 * Convolution of the semiclassical RR cross section with Stobbe corrections
 *  with a Maxwellian characterized by the electron temperature kT yielding the
 *  RR rate coefficient in a plasma as a function of kT. 
*/
void calcAlphaRRplasma(void);



void calc_sigmaRR(void);

void compare_sigmarr(void);

const double Pi = 3.1415926535;  ///< Pi

const double Alpha = 137.0359895;  ///< 1/(fine structure constant) 

const double A0 = 5.29177249e-9; ///< Bohr radius in cm

const double Rydberg=13.6056981; ///< Rydberg's constant in eV

const double Melectron=510998.9;  ///< electron rest mass in eV

const double Matomic = 931494095.4;  ///< atomic mass unit rest mass in eV

const double Clight = 2.99792458E10; ///< speed of Light in cm/s




