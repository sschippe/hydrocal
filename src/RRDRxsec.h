/**
 * @file RRDRxsec.h
 *
 * @brief hydrogenic RR cross sections and DR peak cross sections
 *
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
 * @brief nonrelativistic (dipole approximation) n,l-specific RR cross section (NRD) at zero energy
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrnrd(double z, int n, int l);

///////////////////////////////////////////////////////////////
/**
 * @brief nonrelativistic (dipole approximation) n,l-specific RR cross section (NRD)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrnrd(double eV, double z, int n, int l);

///////////////////////////////////////////////////////////////
/**
 * @brief nonrelativistic (dipole approximation) total RR cross section (NRD of type SIGMARR)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrnrd(double eV, double z, std::vector<double> fraction,
                 int use_fraction, int nmin, int nmax, int lmin);

///////////////////////////////////////////////////////////////
/**
 * @brief relativistic (Dirac, dipole approximation) n,l-specific RR cross section (DID) 
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrdid(double eV, double z, int n, int l);

///////////////////////////////////////////////////////////////
/**
 * @brief relativistic (Dirac, dipole approximation) total RR cross section (DID of type SIGMARR)
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrdid(double eV, double z, std::vector<double> fraction,
                 int use_fraction, int nmin, int nmax, int lmin);

///////////////////////////////////////////////////////////////
/**
 * @brief relativistic (Dirac, full retardation) n,l-specfic RR cross section (DIR) 
 *
 * @return cross section times energy (cm2 eV)
 */
double sigmarrdir(double eV, double z, int n, int l);

///////////////////////////////////////////////////////////////
/**
 * @brief relativistic (Dirac, full retardation) total RR cross section (DIR of type SIGMARR)
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
