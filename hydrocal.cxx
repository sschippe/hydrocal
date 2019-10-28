/** 
    @file hydrocal.cxx
    @brief Main program
 
    @par CREATION  
    @author Stefan Schippers
    @date 1997, 2006, 2007, 2016, 2018
    
    @par VERSION
    @verbatim
    $Id: hydrocal.cxx 631 2019-06-16 06:16:24Z iamp $
    @endverbatim   

   Compilation of hydrogenic atomic structure calculations\n
   and other useful stuff (see hydrocal.cxx for details)\n
                                                             
   Stefan.Schippers@physik.uni-giessen.de\n                                                       
   Atomic and MolecularPhysics\n
   Institute of Experimental Physics I\n
   Justus-Liebig-Universitaet Giessen\n
   Leihgesterner Weg 217, D 35392 Giessen\n
   http://www.uni-giessen.de/amp
 */ 


#include <cstdio>
#include <iostream>
#include "clebsch.h"
#include "lifetime.h"
#include "polari.h"
#include "osci.h"
#include "ratecoef.h"
#include "sigmarr.h"
#include "fraction.h"
#include "rydserie.h"
#include "drmodel.h"
#include "convolute.h"
#include "fieldion.h"
#include "fieldion2.h"
#include "RRphotons.h"
#include "hydromath.h"
#include "kinema.h"
#include "Stark.h"
#include "BetheBloch.h"
#include "autoso1.h"
#include "peakfunctions.h"
#include "stripping.h"
#include "TOPbase.h"
#include "burgess.h"
#include "CoolerMonteCarlo.h"
#include "PIPE-MonteCarlo.h"
#include "svnrevision.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////
/** 
 * mathematical functions (menu item 1)\n
 * <ol>
 * <li>some basic mathematical functions (calls #TestMath)</li>
 * <li>Clebsch Gordan coefficients (calls #testCG)</li>
 * <li>3J symbol (calls #test3J)</li>
 * <li>6J symbol (calls #test6J)</li>
 * <li>9J symbol (calls #test9J)</li>
 * <li>Clebsch Gordan coefficients for transformation to 
 * Stark states (calls #StarkCG)</li>
 * </ol>
 */

void math_group(void)
  {
  int i=1,choice=0;
  char answer;
  do
     {
     printf("\n\n Mathematical functions group: \n\n"); 
     printf(" 1) basic mathematical functions\n");
     printf(" 2) Clebsch Gordan coefficient\n");
     printf(" 3) 3-J Symbol\n");
     printf(" 4) 6-J Symbol\n");
     printf(" 5) 9-J Symbol\n");
     printf(" 6) Clebsch Gordan coefficients for ");
     printf("transformation to Stark states\n");
     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0:  i=0; break;
       case 1:  printf("\n");
                 TestMath();
                 break;
       case 2: printf("\n"); 
	         testCG(); 
                 break;
       case 3: printf("\n");
                 test3J(); 
                 break;
       case 4: printf("\n");
                 test6J(); 
                 break;
       case 5: printf("\n");
                 test9J(); 
                 break;
       case 6: printf("\n");
                 StarkCG(); 
                 break;
      default  : goto error;
      }
    }
  while (i); 
  }

/////////////////////////////////////////////////////////////////////////////
/** 
 * hydrogenic dipole transitions (menu item 2)\n
 * <ol>
 * <li>bound-bound oscillator strengths (calls #testOsciBB)</li>
 * <li>bound-continuum oscillator strengths (calls #testOsciBC)</li>
 * <li>dipole transition rates and lifetimes (calls #testLifetime)</li>
 * <li>relative line strengths (calls #testStrength)</li>
 * </ol>
 */

void dipole_group(void)
  {
  int i=1,choice=0;
  char answer;
  do
     {
     printf("\n\n Hydrogenic dipole transition group: \n\n"); 
     printf(" 1) bound-bound oscillator strength\n");
     printf(" 2) bound-continuum oscillator strength\n");
     printf(" 3) dipole transition rates and radiative lifetimes\n");
     printf(" 4) relative line strength\n"); 
 
     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0:  i=0; break;
       case 1: printf("\n");
                 testOsciBB(); 
                 break;
       case 2: printf("\n");
                 testOsciBC(); 
                 break;
       case 3: printf("\n");
                 testLifetime(); 
                 break;
       case 4: printf("\n");
                 testStrength(); 
                 break;
      default  : goto error;
      }
    }
  while (i); 
  }

/////////////////////////////////////////////////////////////////////////////
/** 
 * (de)polarization calculations (menu item 3)\n
 * <ol>
 * <li> hyperfine slittings (calls #testHF)</li>
 * <li> depolarization factor (calls #testDepol)</li>
 * <li> polarization coefficients (calls #testPol)</li>
 * </ol>
 */

void polari_group(void)
  {
  int i=1,choice=0;
  char answer;
  do
     {
     printf("\n\n Polarization group: \n\n"); 
     printf(" 1) hyperfine splittings\n");
     printf(" 2) depolarization factor\n");
     printf(" 3) polarization coeffcients\n");

     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0:  i=0; break;
       case 1: printf("\n");
                 testHF(); 
                 break;
       case 2: printf("\n");
                 testDepol(); 
                 break;
       case 3: printf("\n");
                 testPol(); 
                 break;
     default  : goto error;
      }
    }
  while (i); 
  }

/////////////////////////////////////////////////////////////////////////////
/** 
 * field ionization calculations (menu item 4)\n
 * <ol>
 * <li> n_gamma by A. Wolf (calls #ngamma)</li>
 * <li> field ionization rates by Damburg \& Kolosov (calls #testFI)</li>
 * <li> surviving nl-fractions (calls #fractionTable)</li>
 * <li> QM field ionization rates (under construction, calls #FIhydrogen)</li>
 * <li> setting up a batch file for running option 3 (calls #setup_batch)</li>
 * <li> delayed cut off (method of Zong et al, calls #fractionTable)</li>
 * </ol>
 */

void fieldion_group(void)
  {
  int i=1,choice=0;
  char answer;
  do
     {
     printf("\n\n Field ionization group: \n\n"); 
     printf(" 1) higest n contributing to RR (ngamma)\n");
     printf(" 2) field ionisation rates and survival probabilities");
     printf(" (Damburg und Kolosov)\n");
     printf(" 3) nl-specific detection probabilities\n");
     printf(" 4) QM field ionisation rates (under construction)\n");
     printf(" 5) setting up batch file for running option 3)\n");
     printf(" 6) nl-specific detection probabilities (method of Zong et al.)\n");
     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0:  i=0; break;
       case 1: printf("\n");
                 ngamma(); 
                 break;
       case 2: printf("\n");
                 testFI(); 
                 break;
       case 3: printf("\n");
                 fractionTable(0); 
                 break;
       case 4: printf("\n");
                 FIhydrogen();
                 break;
       case 5: printf("\n");
                 setup_batch();
                 break;
       case 6: printf("\n");
                 fractionTable(1);
                 break;
     default  : goto error;
      }
    }
  while (i); 
  }

/////////////////////////////////////////////////////////////////////////////
/** 
 *hydrogenic RR cross sections(menu item 5)\n
 * <ol>
 * <li> comparison semi-classical vs. quantum mechanical 
 * (calls #compare_sigmarr)</li>
 * <li> semi-classical or quantum mechanical (calls #calc_sigmaRR)</li>
 * <li> Stobbe correction factors at E=0 (calls #calc_stobbe)</li>
 * <li> Cross sections for photon emission due to RR (calls #calc_emission)</li>
 * </ol>
 */

void sigmarr_group(void)
  {
  int i=1,choice=0;
  char answer;
  do
     {
     printf("\n\n RR cross sections group: \n\n"); 
     printf(" 1) RR cross sections (SCL/QM comparison)\n");
     printf(" 2) RR cross sections (QM or SCL)\n");
     printf(" 3) Stobbe correction factors\n");
     printf(" 4) Photon emssion subsequent to RR \n");

     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0:  i=0; break;
       case 1: printf("\n");
                 compare_sigmarr(); 
                 break;
       case 2: printf("\n");
                 calc_sigmaRR(); 
                 break;
       case 3: printf("\n");
                 calc_stobbe();
                 break;
       case 4: printf("\n");
                 calc_emission();
                 break;
     default  : goto error;
      }
    }
  while (i); 
  }

/////////////////////////////////////////////////////////////////////////////
/** 
 * merged-beams and plasma recombination rate coeffcients (menu item 6)\n
 * <ol>
 * <li> merged beams rate coefficients from RR and DR cross sections
 * (calls #calc_alpha)</li>
 * <li> Monte-Carlo simulation of merged-beams recombination rate coefficients
 * (calls #cooler_MonteCarlo)</li> 
 * <li> plasma rate coefficients from RR and DR cross sections
 * (calls #calc_alpha)</li>
 * <li> plasma rate coefficients from experimental data 
 * (calls #convolute_alpha)</li>
 * </ol>
 */

void ratecoeff_group(void)
  {
  int i=1,choice=0;
  char answer;
  do
     {
     printf("\n\n RR and DR rate coefficients group: \n\n"); 
     printf(" 1) merged-beams rate coefficients from RR and DR cross sections\n");
     printf(" 2) Monte-Carlo simulation of RR and DR merged-beams rate coefficients\n");
     printf(" 3) plasma rate coefficients from RR and DR cross sections\n");
     printf(" 4) plasma rate coefficients from experimental data\n");

     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0:  i=0; break;
       case 1: printf("\n");
	 calc_alpha(1); // ratecoeff.cxx
                 break;
       case 2: printf("\n");
	         cooler_MonteCarlo();  // CoolerMonteCarlo.cxx
                 break;
       case 3: printf("\n");
	 calc_alpha(2);  // ratecoeff.cxx
                 break;
       case 4: printf("\n");
                 convolute_PlasmaRateCoeff();
                 break;
     default  : goto error;
      }
    }
  while (i); 
  }

/////////////////////////////////////////////////////////////////////////////
/** 
 * characteristics of DR spectra (menu item 7)\n
 * <ol>
 * <li>resonance energies of a Rydberg series (calls #rydbergseries)</li>
 * <li>model DR spectrum (calls #calcdr)</li>
 * <li>min E-field for Stark mixing (not implemented, calls #CalcminField)</li>
 * <li>hyperfine splitting of DR resonance (calls #HyperfineDR)</li>
 * </ol>
 */

void modelspe_group(void)
  {
  int i=1,choice=0;
  char answer;
  do
     {
     printf("\n\n DR Model spectrum group: \n\n"); 
     printf(" 1) Resonance energies of a Rydberg series\n");
     printf(" 2) DR model spectrum\n");
     printf(" 3) min E-field for Stark-mixing\n");
     printf(" 4) HF-split DR resonances\n");

     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0:  i=0; break;
       case 1: printf("\n");
                 rydbergseries(); 
                 break;
       case 2: printf("\n");
                 calcdr(); 
                 break;
       case 3: printf("\n");
                // CalcminField(); 
                 break;
       case 4: printf("\n");
                 HyperfineDR(); 
                 break;
      default  : goto error;
      }
    }
  while (i); 
  }

/////////////////////////////////////////////////////////////////////////////
/** 
 * convolutions (menu item 8)\n
 * <ol>
 * <li>merged-beams rate coefficients from RR and DR cross sections
 * (calls #calc_alpha)</li>
 * <li>plasma rate coefficients from RR and DR cross sections
 * (calls #calc_alpha)</li>
 * <li>plasma rate coefficients from experimental data 
 * (calls #convolute_alpha)</li>
 * <li>convolutions of energy dependent cross sections
 * (x,y-data, calls #convflat)</li>
 * <li>convolutions of energy dependent rate coefficients
 * (x,y-data, calls #convflat)</li>
 * <li>sum of peaks (peak parameters, calls #sum_of_peaks)</li>
 * </ol>
 */ 
void convolution_group(void)
  {
  int i=1,choice=0;
  char answer;
  do
     {
     printf("\n\n convolution group: \n\n"); 
     printf(" 1) merged-beams rate coeff. from RR and DR cross sections\n");
     printf(" 2) plasma rate coefficients from RR and DR cross sections\n");
     printf(" 3) plasma rate coefficients from experimental data\n");
     printf(" 4) convolutions of energy dependent cross sections\n");
     printf("    contained as x-y columns in a file\n");
     printf(" 5) convolutions of energy dependent rate coeff.\n");
     printf("    contained as x-y columns in a file\n");
     printf(" 6) sum of peaks with peak parameters from a file\n");
     printf(" 7) photoionization cross sections from TOPbase\n");

     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0:  i=0; break;
       case 1: printf("\n");
	 calc_alpha(1); // ratecoeff.cxx
                 break;
       case 2: printf("\n");
	 calc_alpha(2); // ratecoeff.cxx
                 break;
       case 3: printf("\n");
	 convolute_PlasmaRateCoeff(); // convolute.cxx
                 break;
       case 4: printf("\n");
	 convolute_Xsec_MBrateCoeff(1);  //convolute.cxx
                 break;
       case 5: printf("\n");
	 convolute_Xsec_MBrateCoeff(2); // convolute.cxx
                 break;
       case 6: printf("\n");
	 sum_of_peaks(); // peakfunctions.cxx
                 break;
       case 7: printf("\n");
	 convolute_TOPbase_PI(); // TOPbase.cxx
                 break;
     default  : goto error;
      }
    }
  while (i); 
  }

/////////////////////////////////////////////////////////////////////////////
/** 
 * Monte Carlo simulations (menu item 9)\n
 * <ol>
 * <li>Molecular breakup in the PIPE setup
 * (calls #PIPE-MonteCarlo)</li>
 * <li> Monte-Carlo simulation of merged-beams recombination rate coefficients
 * (calls #cooler_MonteCarlo)</li> 
 * </ol>
 */ 
void Monte_Carlo_group(void)
  {
  int i=1,choice=0;
  char answer;
  do
     {
     printf("\n\n Monte Carlo group: \n\n"); 
     printf(" 1) Monte-Carlo simulation of molecular breakup in the PIPE setup\n");
     printf(" 2) Monte-Carlo simulation of RR and DR merged-beams rate coefficients\n");

     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0:  i=0; break;
       case 1: printf("\n");
#ifdef EIGEN
	         PIPE_MonteCarlo();  // PIPE-MonteCarlo.cxx
#else
                 cout << endl << "*******************************************************************";
		 cout << endl << "*** This option reqires the package libeigen3-dev.              ***";
                 cout << endl << "*** If available recompile hydrogen using 'make hydrocal-eigen' ***";
                 cout << endl << "*******************************************************************";
		 cout << endl;  
#endif
                 break;
       case 2: printf("\n");
	         cooler_MonteCarlo();  // CoolerMonteCarlo.cxx
                 break;
     default  : goto error;
      }
    }
  while (i); 
  }

/////////////////////////////////////////////////////////////////////////////
/** 
 * miscellaneous calculations (menu item 10)\n
 * <ol>
 * <li>relativistic transformation for merged beams (calls #kinema)</li>
 * <li>Weizsäcker nuclear mass formula (calls #TestWeizmass)</li>
 * <li>temperature average of v1 in rate enhancement 
 * (experimental, calls #calc_v1_averaged)</li>
 * <li>Bethe-Bloch stooping power (calls #BetheBloch)</li>
 * <li>reading and processing of autostructure data (calls #TestAutosO1)</li>
 * <li>Burgess formula for PI cross sections (calls #BurgessPIxsec)</li>
 * </ol>
 */ 
void miscellaneus_group(void)
  {
  int i=1,choice=0;
  char answer;
  do
     {
     printf("\n\n calculate this and that: \n\n"); 
     printf(" 1) Relativistic transformation for merged beams\n");
     printf(" 2) Weizsaecker mass formula\n");
     printf(" 3) Bethe Bloch\n");
     printf(" 4) autostructre o1 read\n");
     printf(" 5) charge state distribution after stripping\n");
     printf(" 6) PI cross sections (Burgess formula)\n");

     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0:  i=0; break;
       case 1: printf("\n");
                 kinema(); 
                 break;
       case 2: printf("\n");
                 TestWeizmass(); 
                 break;
       case 3: printf("\n");
                 BetheBloch();
	         break;
       case 4: printf("\n");
	         TestAutosO1();
	         break;
       case 5: printf("\n");
	         stripping();
	         break;
       case 6: printf("\n");
	         BurgessPIxsec();
	         break;
      default  : goto error;
      }
     choice=0;
    }
  while (i); 
  }


/////////////////////////////////////////////////////////////////////////////
/**
 * the main program from where the above menues can be accessed
 */
int  main(void)
  {
  int i=1,choice=0;
  char answer;
  printf(" **************************************************************\n");
  printf(" *  hydrocal, revision %-20s                   *\n",SVNrevision);
  printf(" *                                                            *\n");
  printf(" *         Compilation of hydrogenic atomic structure         *\n");
  printf(" *               and other useful calculations                *\n");
  printf(" *                                                            *\n");
  printf(" *               Stefan Schippers, 1997...2019                *\n");
  printf(" *                                                            *\n");
  printf(" *            Justus-Liebig-Universitaet Giessen              *\n"); 
  printf(" *            I. Physikalisches Institut                      *\n");
  printf(" *            Atom- und Molekuelphysik                        *\n");
  printf(" *            Leihgesterner Weg 217, D 35392 Giessen          *\n");
  printf(" *                                                            *\n");
  printf(" *        e-mail: Stefan.Schippers@physik.uni-giessen.de      *\n");
  printf(" **************************************************************\n");

  do
    {
     printf("\n\n Calculate what ? \n\n"); 
     printf(" 1) mathematical functions and 3n-j symbols\n");
     printf(" 2) hydrogenic dipole transitions\n");
     printf(" 3) hyperfine splittings and polarisations\n");
     printf(" 4) field ionization\n");
     printf(" 5) RR cross sections\n");
     printf(" 6) RR and DR rate coefficients\n");
     printf(" 7) model DR spectra\n");
     printf(" 8) convolutions\n");
     printf(" 9) Monte Carlo simulations\n");
     printf("10) this and that\n");
     error: printf("\n Make your choice (0 quits )               : ");
     scanf("%d",&choice); 
     switch (choice)
      {
       case 0: i=0;
                 break;
       case 1: printf("\n");
                 math_group(); 
                 break;
       case 2: printf("\n"); 
                 dipole_group(); 
                 break;
       case 3: printf("\n");
                 polari_group(); 
                 break;
       case 4: printf("\n");
                 fieldion_group(); 
                 break;
       case 5: printf("\n");
                 sigmarr_group(); 
                 break;
       case 6: printf("\n");
                 ratecoeff_group(); 
                 break;
       case 7: printf("\n");
                 modelspe_group(); 
                 break;
       case 8: printf("\n");
                 convolution_group();
                 break;
       case 9: printf("\n");
                 Monte_Carlo_group();
                 break;
       case 10: printf("\n");
                 miscellaneus_group();
                 break;
      default  : goto error;
      }
      choice=0;
    }
  while (i); 
  printf("\n");
  printf("\n");
 return 1;
}







