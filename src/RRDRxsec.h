/**
 * @file RRDRxsec.h
 *
 * @brief hydrogenic RR cross sections and DR peak cross sections
 *
 * $Id: RRDRxsec.h 2039 2026-07-20 07:57:32Z iamp $
// SPDX-License-Identifier: MIT
 */

#pragma once

#include "fele.h"
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief type definition allowing for passing cross section functions as parameters
 */
typedef double (*SIGMARR)(double /* e */, double /* z */, std::vector<double> /* faction */, int /* use_faction */,
			  int /* nmin */ , int /*nmax */ , int /* lmin */);


///////////////////////////////////////////////////////////////
/**
 * @brief quantum mechanical PI cross section
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmaPIqm(double Eph_eV, double z, int nmin, int n, int l, double IP);

///////////////////////////////////////////////////////////////
/**
 *  @brief quantum mechanical PI cross section (of type SIGMARR)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmaPIqm(double eV, double z, std::vector<double> fraction,
                 int use_fraction, int nmin, int nmax, int lmin);

///////////////////////////////////////////////////////////////
/**
 * @brief quantum mechanical RR cross section times energy for energy->0
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrqm(double z, int n, int l);

///////////////////////////////////////////////////////////////
/**
 * @brief quantum mechanical RR cross section
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrqm(double eV, double z, int n, int l);

///////////////////////////////////////////////////////////////
/**
 * @brief quantum mechanical RR cross section (of type SIGMARR)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrqm(double eV, double z, std::vector<double> fraction,
                 int use_fraction, int nmin, int nmax, int lmin);

///////////////////////////////////////////////////////////////
/**
 * @brief relativistic (Dirac) RR cross section 
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrdir(double eV, double z, int n, int l);

///////////////////////////////////////////////////////////////
/**
 * @brief relativistic (Dirac) RR cross section (of type SIGMARR)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrdir(double eV, double z, std::vector<double> fraction,
                 int use_fraction, int nmin, int nmax, int lmin);

///////////////////////////////////////////////////////////////
/**
 * @brief Stobbe correction factor for semiclassical RR cross section
 */
double stobbe(int n);

///////////////////////////////////////////////////////////////
/**
 * @brief calculation of Stobbe correction factors
 */
void calc_stobbe(void);

///////////////////////////////////////////////////////////////
/**
 * @returns Stobbe correction factor for specified atomic shell
 */
double kstobbe(int n);

///////////////////////////////////////////////////////////////
/**
 * @brief semiclassical RR cross section (of type SIGMARR)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrscl(double eV, double z, std::vector<double> fraction,
                  int use_fraction, int nmin, int nmax, int lmin);

///////////////////////////////////////////////////////////////
/**
 * @brief semiclassical RR cross section with high-n integration (of type SIGMARR)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrscl2(double eV, double z, std::vector<double> fraction,
                   int use_fraction, int nmin, int nmax, int lmin);

///////////////////////////////////////////////////////////////
/**
 * @brief quantum mechanical RR cross section averaged over angular momentum  (of type SIGMARR)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarraqm(double eV, double z, std::vector<double> fraction,
                  int use_fraction, int nmin, int nmax, int lmin);

///////////////////////////////////////////////////////////////
/**
 * @brief test cross section for comparisons (of type SIGMARR)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrtest(double e, double z, std::vector<double> fraction,
                   int use_fraction, int nmin, int nmax, int lmin);


///////////////////////////////////////////////////////////////
/**
 * @brief convolution of a lorentzian resonance cross section with electron energy distribution
 * 
 * @return cross section times energy (cm2 eV)
 */
double sigmalorentzian(double e, double er, std::vector<double> fraction,
                       int use_fraction, int nmin, int nmax, int lmin);


///////////////////////////////////////////////////////////////
/**
 * @brief convolution of a delta peak resonance cross section with electron energy distribution
 * 
 * @return merged-beams rate coefficient
 */
double deltapeak(FELE fele, double er, double eV, double ktpar, double ktperp, double a);

///////////////////////////////////////////////////////////////
/**
 * @brief plasma rate coefficient for RR from semiclassical cross sections
 */
double alphaRRplasmaSCL(double kt, double z, int nmin, int nmax, int nele);

///////////////////////////////////////////////////////////////
/**
 * @brief plasma cooling coefficient for RR from semiclassical cross sections
 */
double betaRRplasmaSCL(double kt, double z, int nmin, int nmax, int nele);

///////////////////////////////////////////////////////////////
/**
 * @brief calculation of interstellar (plasma) RR rate coefficients
 */
void calc_sigmaRR(void);

///////////////////////////////////////////////////////////////
/**
 * @brief calculate RR rate coefficient from cross section stored in file
 */
void calcAlphaRRplasma(void);

///////////////////////////////////////////////////////////////
/** 
 * @brief initializes arrays for interpolation
 */
void init_sigma_interpolation(std::vector<double> energy,
                              std::vector<double> xsec, int npts);
///////////////////////////////////////////////////////////////
/**
 * @brief interpolated cross section(of type SIGMARR)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigma_interpolated(double e, double z, std::vector<double> fraction,
                          int use_fraction, int nmin, int nmax, int lmin);
