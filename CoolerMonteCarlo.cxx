/**
 * @file CoolerMonteCarlo.cxx
 *
 * @brief Monte-Carlo simulation of merged-beams recombination rate coefficients
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: CoolerMonteCarlo.cxx 614 2019-03-29 14:01:00Z iamp $
   @endverbatim
 *
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <chrono>     // clocks and time
#include <random>     // random number generators and distributions  
#include <functional> // for "bind (random_generator, random_distribution)"
#include "readxsec.h"
#include "sigmarr.h"
#include "svnrevision.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief relativistic beta from kinetic energy
 *
 * @param Ekin kinetic energy
 *
 * @return beta velocity divided by speed of light
 */
double RelBeta(double Ekin)
{
  double gamma = 1.0+Ekin/Melectron;
  return sqrt(1.0-1.0/(gamma*gamma));
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculates the electron-ion collision energy from the LAB electron energy the and the LAB ion energy
 *
 * @param Erel electron-ion collision energy in eV
 * @param Ecool cooling energy in lab frame in eV
 * @param Aion ion mass number
 * @param costheta cosine of angle between electron beam and ion beam
 *
 * @return electron-ion collision energy in center-of-mass frame
 */
double LABtoCM(double Eele, double Ecool, double Aion, double costheta=1.0)
{
  double mass_ratio, mass_factor, gamma_ele, gamma_ion, Gamma;
  // Melectron and Matomic defined in sigmarr.h
  mass_ratio = Melectron/(Aion*Matomic);  // Melectron and Matomic defined in sigmarr.h
  mass_factor = 2.0*mass_ratio/((1.0+mass_ratio)*(1.0+mass_ratio));
  gamma_ele = 1.0+Eele/Melectron;
  gamma_ion = 1.0+Ecool/Melectron;
  Gamma = gamma_ion*gamma_ele-sqrt((gamma_ion*gamma_ion-1.0)*(gamma_ele*gamma_ele-1.0))*costheta;
  //cout << gamma_ion << " " << gamma_ele << " " << Gamma-1.0 << " " << sqrt(1.0+mass_factor*(Gamma-1.0))-1.0 << endl;
  return (Aion*Matomic+Melectron)*(sqrt(1.0+mass_factor*(Gamma-1.0))-1.0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculates the electron energy from the electron-ion collision energy and the ion energy
 *
 * Approximate formula neglecting Eree*Erel/(mec2*mic2) term in transformation, valid for Erel << 30 MeV sqrt(A)
 * see appendic of PhD thesis of Wolgang Spies (JLU Giessen, 1995)
 *
 * @param Erel electron-ion collision energy in eV
 * @param Ecool cooling energy in lab frame in eV
 * @param Aion ion mass number
 * @param Eele1 on exit: first solution for electron energy in lab frame
 * @param Eele2 on exit: second solution for electron energy in lab frame
 */
void CMtoLABele(double Erel, double Ecool, double Aion, double &Eele1, double &Eele2)
{
  // Melectron and Matomic defined in sigmarr.h
  double Mreduced, gamma_ele, gamma_ion, gamma_rel, sqrt_term;
  Mreduced  = Melectron*Aion*Matomic/(Melectron+Aion*Matomic);
  gamma_rel = 1.0+Erel/Mreduced;
  gamma_ion = 1.0+Ecool/Melectron;
  sqrt_term = sqrt((gamma_rel*gamma_rel-1.0)*(gamma_ion*gamma_ion-1.0));
  gamma_ele = gamma_rel*gamma_ion-sqrt_term;
  Eele1 = Melectron*(gamma_ele-1.0);
  gamma_ele = gamma_rel*gamma_ion+sqrt_term;
  Eele2 = Melectron*(gamma_ele-1.0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
bool CoolerAngleESR(double zpos, double &angle); // coded at the end of this file 

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *  @brief calculation of the 3D electron velocity in the drift tube 
 *
 *  Takes fringe potentials at the entrance and at the exit of the drift tube into account,
 *  as well as angles due to the magnetic guiding field in the cooler
 *
 * @param z position in the cooler ( 0.0 < z < 250.0 )
 * @param Ucath cathode voltage
 * @param eLAB electron energy in the laboratory frame
 * @param mode 0: flat potential, no angles; 1: drift tube potential, no angles; 2: drift-tube potential with angles
 * @param offset_angle overall offset angle in mrad subtracted at all positions
 * @param ele_beta_x on exit: x-component of electron velocity at position z
 * @param ele_beta_y on exit: y-component of electron velocity at position z
 * @param ele_beta_z on exit: z-component of electron velocity at position z
 *
 * @return false, if z is out of range
 * @return true, if z is ok
 */ 
bool DriftTube(double z, double Ucath, double eLAB, int mode, double offset_angle, double &ele_beta_x, double &ele_beta_y, double &ele_beta_z)
{
  const double cooler_length = 250.0; // in cm
  const double drifttube_halflength = 97.0;  // in cm
  const double drifttube_radius =  10.0;  // in cm

  double theta; 
  bool error_flag = CoolerAngleESR(z,theta); // returns theta in rad
  if (error_flag ) return true; // z is out of range
    
  double Udrift=eLAB-Ucath;
  double ele_beta;
  if (mode==0)
    { // flat potential distribution, i.e., not fringe fields
      ele_beta= RelBeta(Ucath+Udrift);
    }
  else
    { // drift tube potential with fringe fields
      double zz = z-0.5*cooler_length; // zz=0 in the center of the cooler
      if (zz <= 0.0)
	{
	  Udrift *= 0.5+0.5*tanh(1.318*(zz+drifttube_halflength)/drifttube_radius);
	}
      else
	{
	  Udrift *= 0.5-0.5*tanh(1.318*(zz-drifttube_halflength)/drifttube_radius);
	}
	ele_beta= RelBeta(Ucath+Udrift);
     }  

  if (mode==2)
    { // take angles into account
      theta -= offset_angle*0.001; // in rad
      ele_beta_z = ele_beta*cos(theta);   
      ele_beta_y = ele_beta*sin(theta);
    }
  else
    {
      ele_beta_z = ele_beta;   
      ele_beta_y = 0.0; 
    }
  ele_beta_x = 0.0; 
  return false; // normal exit
}

//////////////////////////////////////////////////////////////////
/**
 * @brief writes drift-tube potentials and angles as function of position to a file
 *
 * for inspection only, nomally not used
 */
void PrintDriftTubeData()
{
  const double cooler_length = 250.0; // in cm
  const double drifttube_halflength = 97.0;  // in cm
  const double drifttube_radius =  10.0;  // in cm

  fstream fout;
  fout.open("drifttube.dat",fstream::out);
  for (int i=0; i<=2500; i++)
    {
      double z= 0.1*(double)i;
      double theta; 
      bool error_flag = CoolerAngleESR(z,theta); // returns theta in rad
      if (error_flag) continue;
      double zz = z-0.5*cooler_length; // zz=0 in the center of the cooler
      double Udrift=1.0;  
      if (zz <= 0.0)
	{
	  Udrift *= 0.5+0.5*tanh(1.318*(zz+drifttube_halflength)/drifttube_radius);
	}
      else
	{
	  Udrift *= 0.5-0.5*tanh(1.318*(zz-drifttube_halflength)/drifttube_radius);
	}
      fout << z << ", " << Udrift << ", " << theta << endl;
    }
  fout.close();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief binning of Lorentzian DR cross sections
 *
 * Analytically integrates DR cross section over histogram bins
 *
 * Note that inconsistencies can arrise if the DR cross section
 * is already binned (e.g. by AUTOSTRUCTRUE) and the original 
 * bin size is larger than the one used here. 
 *
 * @param nDR number of DR resonances
 * @param nbin number of histogram bins 
 * @param bin_width uniform width of all bins in eV
 * @param eHmin min energy of histogram in eV
 * @param eDR array of DR resonance energies in eV
 * @param sDR array of DR resonance strengths in cm2 eV  
 * @param wDR array of DR resonance widths in eV
 * @param sigmaDR on exit binned cross section (averaged per bin) in cm2  
 *
 * @return sum of binned resonance strengths
 */
double BinSigmaDR(int nDR, int nbin, double bin_width, double eHmin, 
		double *eDR, double *sDR, double *wDR, double *sigmaDR)
{
  const double eps = 1E-9;
  int n,i;
  for (i=0; i<nbin; i++) { sigmaDR[i] = 0.0; }
  for (i=0; i<nbin; i++)
    {
      double emin = eHmin+i*bin_width; // min energy of bin
      double emax = emin+bin_width;     // max energy of bin
      for (n=0; n<nDR; n++)
	{
	  double eres = eDR[n];
	  double wres = fabs(wDR[n]);
	  if ( wres < eps )  
	    { // treat narrow resonances as delta-like resonances
	      if  ( (eres >= emin) && (eres < emax) ) sigmaDR[i] += sDR[n];
	    }
	  else
	    { // integrate the Lorentzian resonance from lower to upper bin boundary
	      double xmin = 2.0*(emin-eres)/wres;
	      double xmax = 2.0*(emax-eres)/wres;
	      sigmaDR[i] += sDR[n]*(atan(xmax)-atan(xmin))/Pi;
	    }
	}//end for(n...
    } // end for(i...
  double sumbins=0.0;
  for (i=0; i<nbin; i++) { sumbins += sigmaDR[i]; sigmaDR[i]/=bin_width; }
  return sumbins;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Entry to Monte Carlo simulation of electron-ion merged-beams rate coefficients
 */
int cooler_MonteCarlo()
{
  cout<< endl << " Monte Carlo simulation of RR and DR spectra from electron-ion merged-beams arrangements in ion-storage rings";
  cout<< endl;
  //PrintDriftTubeData(); // for inspection only
  
  int nDR = 1, theo_mode = 0, number_of_levels = 1, level_number = 1;
  char theo_filename[200];   

  info_DRtheory();
  while ((theo_mode==0)||(theo_mode>11))
    {
      cout << endl << " Give name of DR theory data file ......................: ";
      cin >> theo_filename;
      nDR = open_DRtheory(theo_filename, theo_mode, number_of_levels);
      if (nDR<0)
	{
	  cout << endl << " File " << theo_filename << " not found." << endl;
	  theo_mode=0;
	}
    }
  cout << endl << " " << nDR << " DR resonances found." << endl; 
  if (nDR==0) return 1;
  if (number_of_levels>1)
    {
      cout << endl << " Data for " << number_of_levels << " different initial levels found.";
      cout << endl << " Give the number (in the range 1 - "<<number_of_levels<<") of the level to be used: ";
      cin >> level_number;
      if ( (level_number<1) || (level_number>number_of_levels) )
	{
	  cout << endl << "\n ERROR: Level number out of range!" << endl << endl;
	  exit(0);
	}
    }
  double *eDR = new double[nDR];
  double *sDR = new double[nDR];
  double *wDR = new double[nDR];
  double eDRmin, eDRmax;
  read_DRtheory(theo_filename,theo_mode,level_number,eDR,sDR,wDR,nDR,eDRmin,eDRmax);
  cout << endl << " DR resonances ranging from " << eDRmin << " eV to " << eDRmax << " eV" << endl;
  double totalDRstrength=0.0;
  for (int n=0; n<nDR; n++) { totalDRstrength += sDR[n];}

  int use_fraction=0, nint=0, RRnmin, RRnmax, RRlmin, RRnele;
  double RRzeff=92.0;
  cout << endl << " Give Zeff, nmin, lmin, nmax for RR.....................: ";
  cin >> RRzeff >> RRnmin >> RRlmin >> RRnmax;
  cout << endl << " Give number of electrons in min n,l-subshell ..........: ";
  cin >> RRnele;
  double RRfraction[2]; // we only need RRfraction[0];
  RRfraction[0] = 1.0 - RRnele/(4.0*RRlmin+2.0);
 
  double ele_ktpar, ele_ktperp;
  cout << endl << " Give electron beam temperatures kTpar and kTperp (meV) : ";
  cin >> ele_ktpar >> ele_ktperp;
  ele_ktpar *= 0.001;   // meV -> eV
  ele_ktperp *= 0.001;  // meV -> eV

  double Ucool, Ucath, ion_mass, ion_deltap, ion_ktperp;
  cout << endl << " Give cooling voltage (V) ..............................: ";
  cin >> Ucool;
  cout << endl << " Give cathode voltage (V) ..............................: ";
  cin >>Ucath;
  cout << endl << " Give ion mass (u) .....................................: ";
  cin >> ion_mass;
  cout << endl << " Give relative momentum deltaP/P (FWHM) of ion beam ....: ";
  cin >> ion_deltap;
  cout << endl << " Give transverse ion temperature (meV) .................: ";
  cin >> ion_ktperp;
  ion_ktperp *= 0.001;  // meV -> eV

  int drift_tube_mode;
  cout << endl << " Choice of drift-tube modes";
  cout << endl << "     mode 0: flat potential, no angles";
  cout << endl << "     mode 1: drift-tube potential, no angles";
  cout << endl << "     mode 2: drift-tube potential with angles";
  cout << endl << " Make a choice (0/1/2) .................................: ";	 
  cin >> drift_tube_mode;
  if (drift_tube_mode<0) drift_tube_mode=0;
  if (drift_tube_mode>2) drift_tube_mode=2;
  double offset_angle = 0.0;
  cout << endl << " Give additional offset angle (mrad) ...................: ";
  cin >> offset_angle;
  if (drift_tube_mode<2) offset_angle = 0.0;

  // **************************** determine relevant energy range ********************************
  // *
  // * We want to step through LAB electron energies with equidistant steps.
  // * To this end, we calculate the LAB electron energy range from the CM energy range of interest
  // *
  double eCMspread = sqrt(ele_ktperp*log(2.0)*ele_ktperp*log(2.0)+16.0*log(2)*ele_ktpar*eDRmin);
  cout << endl << endl << " Minmum electron energy spread in CM system .: " << eCMspread << " eV";   
  double eCMmin   = eDRmin-10.0*eCMspread;   if (eCMmin <= 0.0) { eCMmin = 0.1;}
  double eCMmax   = eDRmax+10.0*eCMspread;
  double eCMdelta = 0.1*eCMspread;
  double eLAB1min, eLAB1max, eLAB1delta,eLAB2min,eLAB2max, eLAB2delta; 
  CMtoLABele(eCMmin,Ucool,ion_mass,eLAB1min,eLAB2min);
  CMtoLABele(eCMmin+eCMdelta,Ucool,ion_mass,eLAB1delta,eLAB2delta);  
  eLAB1delta -= eLAB1min;
  eLAB2delta -= eLAB2min;
  CMtoLABele(eCMmax,Ucool,ion_mass,eLAB1max,eLAB2max);
  if (eLAB1delta<0) { double swap=eLAB1min; eLAB1min=eLAB1max; eLAB1max=swap; eLAB1delta = -eLAB1delta; }
  if (eLAB2delta<0) { double swap=eLAB2min; eLAB2min=eLAB2max; eLAB2max=swap; eLAB2delta = -eLAB2delta; }
  int nLAB1steps = (int)(1.1+(eLAB1max-eLAB1min)/eLAB1delta);
  int nLAB2steps = (int)(1.1+(eLAB2max-eLAB2min)/eLAB2delta);
  cout << endl << " Suggested CM energy range from DR resonances: " << eCMmin << " - " << eCMmax << " eV";
  cout << endl << " Suggested smallest CM energy step size......: " << eCMdelta << " eV";
  cout << endl << " Resulting LAB electron energy ranges  (step sizes) and numbers of points";
  cout << endl << "    " << eLAB1min << " - " << eLAB1max << " (" << eLAB1delta << ")  :" << nLAB1steps;  
  cout << endl << "    " << eLAB2min << " - " << eLAB2max << " (" << eLAB2delta << ")  :" << nLAB2steps;  
  cout << endl;
  double eLABmin, eLABmax,eLABdelta;
  char answer;
  cout << endl << " Full calculation or energy distribution ? (f/d) .......: ";
  cin >> answer; 
  bool distri_flag = ( (answer=='d') || (answer=='D') );
  cout << endl << " Give min LAB electron energy (eV) .....................: ";
  cin >> eLABmin;
  cout << endl << " Give max LAB electron energy (eV) .....................: ";
  cin >> eLABmax;
  cout << endl << " Give LAB electron energy step size (eV)................: ";
  cin >> eLABdelta;
  if (distri_flag)
    {
      cout << endl << " The energy distribution is only calculated for the max LAB energy." << endl;
      eLABmin = eLABmax;
    }
  int nLABsteps = distri_flag ? 1 : (int)(1.1+(eLABmax-eLABmin)/eLABdelta);
  
  // ******************  seed for the  random generator  ***************************
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
 
  double bin_width_factor;
  cout << endl << " Give bin-width factor (normally 1.0) ..................: ";
  cin >> bin_width_factor;
  int MonteCarlo_iterations;
  cout << endl << " Give the number of Monte-Carlo iterations .............: ";
  cin >> MonteCarlo_iterations;
  
  // ********************* open output file and write header ************************         
  fstream fout;
  char outfilename[200], extension[5];
  if (distri_flag) { sprintf(extension,".mcd"); } else { sprintf(extension,".mcr"); }
  cout << endl << " Give filename for output (" << extension << ") .......................: ";
  cin >> outfilename;
  strcat(outfilename,extension);
  
  fout.open(outfilename,fstream::out);
  fout << "#######################################################################################################################" << endl;
  if (distri_flag)
    {
      fout << "### CM-energy distribution from a Monte-Carlo simulation of an electron-ion merged-beams experiment" << endl;
    }
  else
    {
      fout << "### Rate coefficients from a Monte-Carlo simulation of an electron-ion merged-beams experiment" << endl;
    }
  fout << "###" << endl;
  fout << "###  hydrocal SVN revision : " << SVNrevision << endl; // defined in SVNrevision.h
  fout << "###               filename : " << outfilename << endl;
  fout << "###            DR filename : " << theo_filename << endl;
  if (distri_flag)
    {
      fout << "###   lab elec energy (eV) : " << eLABmin << endl;
      fout << "###   relative energy (eV) : " << LABtoCM(eLABmin,Ucool,ion_mass) << endl;
    }
  else
    {
      fout << "###            RR Zeff     : " << RRzeff << endl;
      fout << "###            RR nmin     : " << RRnmin << endl;
      fout << "###            RR lmin     : " << RRlmin << endl;
      fout << "###            RR nele     : " << RRnele << endl;
      fout << "###            RR nmax     : " << RRnmax << endl;
    }
  fout << "###  electron kTpar  (meV) : " << ele_ktpar*1000.0 << endl;
  fout << "###  electron kTperp (meV) : " << ele_ktperp*1000.0 << endl;
  fout << "###   cooling voltage (eV) : " << Ucool << endl;
  fout << "###   cathode voltage (eV) : " << Ucath << endl;
  fout << "###       ion mass (u)     : " << ion_mass << endl;
  fout << "###       ion delta_p/p    : " << ion_deltap << endl;
  fout << "###       ion kTperp (meV) : " << ion_ktperp*1000.0 << endl;
  fout << "###        drift-tube mode : " << drift_tube_mode << endl;
  fout << "###    offset angle (mrad) : " << offset_angle << endl;
  fout << "###            random seed : " << seed << endl;
  fout << "###       bin-width factor : " << bin_width_factor << endl;
  fout << "### Monte-Carlo iterations : " << MonteCarlo_iterations << endl;
  if (distri_flag)
    {
      fout << "###################################################################################################" << endl;
      fout << "### CM energy (eV) [1]     CM-energy distribution (1/eV) [2]      binned DR cross section (cm2) [3]" << endl;
      fout << "###------------------------------------------------------------------------------------------------" << endl;
    }
  else
    {
      fout << "#############################################################################################################################" << endl;
      fout << "### LAB energy (eV) [1]   relative energy (eV) [2]    alpha DR (cm3/s) [3]    alpha RR (cm3/s) [4]    alpha RR+DR (cm3/s) [5]" << endl;
      fout << "###--------------------------------------------------------------------------------------------------------------------------" << endl;
    }

  cout << endl << " Additional diagnostic output on screen? (y/n) .........: ";
  cin >> answer; 
  bool diagnostic_flag = ( (answer=='y') || (answer=='y') );

  cout << endl << endl << " now calculating, be patient " << flush;

  //************************ initialize random generators and random distributions  *******************

  mt19937_64 randgen(seed);  // 64bit "Mersenne Twister 19937" random number generator
  const double ln2=log(2.0);
  double ion_beta = RelBeta(Ucool);

  uniform_real_distribution<double> dist_cooler(2.0,248.0); // overlap length of 246 cm, magnetic field angles are provided for 0-250 cm 
  normal_distribution<double> dist_ele_perp(0.0,sqrt(ele_ktperp/Melectron));
  normal_distribution<double> dist_ele_par(0.0,sqrt(ele_ktpar/Melectron));
  normal_distribution<double> dist_ion_perp(0.0,sqrt(ion_ktperp/(ion_mass*Matomic)));
  normal_distribution<double> dist_ion_par(0.0,0.25*ion_deltap*ion_beta/sqrt(ln2)); // implies conversion from FWHM to Gaussian sigma
  //  Melectron and Matomic defined in sigmarr.h
  double Mreduced  = Melectron*ion_mass*Matomic/(Melectron+ion_mass*Matomic);

  // *************************************  loop over equidistant lab energies  ************************
  for (int n=0; n < nLABsteps; n++)
    {
      double eLAB = eLABmin+n*eLABdelta;
      double erel = LABtoCM(eLAB,Ucool,ion_mass);
      double erel_spread = sqrt(ele_ktperp*ele_ktperp*ln2*ln2+16.0*ln2*ele_ktpar*erel);

      // ******************** preparation of the discrete electron energy distribution *****************
      
      // We want to sort the random CM energies into a histrogram (dist_eCM).
      // The bin width of the histrogram is derived from the expected energy spread.
      double bin_width = 0.1*erel_spread*bin_width_factor;
      double eCMmin = 0, eCMmax = erel+20.0*erel_spread; // energy interval of histrogram
      int nbin = (int)round((eCMmax-eCMmin)/bin_width);  // number of bins
      // The number of bins is not very costly in the convolution below, but may be in the
      // binning of the DR cross sections. A possible future improvement could be to reuse
      // the same bins for a number of LAB energies.
      double *dist_eCM = new double[nbin];
      double *sigmaDR = new double[nbin];
      for (int j=0; j<nbin; j++) { dist_eCM[j] = 0.0; sigmaDR[j] = 0.0; }
      // average the DR cross section over the bins
      double binnedDRstrength = BinSigmaDR(nDR, nbin, bin_width, eCMmin, eDR, sDR, wDR, sigmaDR);
      if (diagnostic_flag)
	{
	  cout << endl << " Elab: " << eLAB << " eV, Erel: " << erel << " eV, "<< nbin << " bins, bin width: " << bin_width;
	  cout << " eV, binned DR strength: " << binnedDRstrength << " eV cm2, total DR strength: ";
	  cout << totalDRstrength << " eV cm2";
	}
      double dist_eCM_count = 0.0;
      
      // **************************  Monte-Carlo simulation of the CM energy distribution **************
      
      for (int i=0; i<MonteCarlo_iterations;i++)
	{ // z-axis in beam direction
          double zpos = dist_cooler(randgen);
          double ele_beta_x=0.0, ele_beta_y=0.0, ele_beta_z=0.0;          
          if (DriftTube(zpos,Ucath,eLAB,drift_tube_mode,offset_angle,ele_beta_x,ele_beta_y,ele_beta_z)) continue; // if zpos out of range
          // add velocity components due to beam temperatures
	  ele_beta_x += dist_ele_perp(randgen);
	  ele_beta_y += dist_ele_perp(randgen);
	  ele_beta_z += dist_ele_par(randgen);
	  double ion_beta_x = dist_ion_perp(randgen);
	  double ion_beta_y = dist_ion_perp(randgen);
	  double ion_beta_z = ion_beta+dist_ion_par(randgen);          	  
	  // Relativistic addition of velocities, we take for granted that the z-component is much larger than the x- and y-components.
	  // A thorough treatment would require matrix transformations (Lorentz boosts and rotations).
	  double gamma = 1.0/sqrt(1-ion_beta_z*ion_beta_z); // approximate gamma
	  double prod_beta = ion_beta_x*ele_beta_x+ion_beta_y*ele_beta_y+ion_beta_z*ele_beta_z; // scalar product 
          double rel_beta_z = (ele_beta_z-ion_beta_z)/(1-prod_beta);
	  double rel_beta_x = (ele_beta_x-ion_beta_x)/(1-ion_beta_z*ele_beta_z)/gamma;
	  double rel_beta_y = (ele_beta_y-ion_beta_y)/(1-ion_beta_z*ele_beta_z)/gamma;
	  double rel_beta_sqr = rel_beta_x*rel_beta_x+rel_beta_y*rel_beta_y+rel_beta_z*rel_beta_z; 
	  gamma = 1.0/sqrt(1.0-rel_beta_sqr); // better gamma
	  double eCM = Mreduced*(gamma-1.0);
	  int bin =  (int)floor((eCM-eCMmin)/bin_width);
	  if ((bin>=0) && (bin<nbin))
	    {
	      dist_eCM[bin] += 1.0;
	      dist_eCM_count += 1.0;
	    }
        } // end for(i ...) end of Monte-Carlo loop
      //*******************************************

      //****************  Convolution of DR and RR cross sections with the CM energy distribution *********

      double alphaDR=0.0, alphaRR=0.0;
      for (int j=0; j<nbin; j++)
	{
	  dist_eCM[j] /= dist_eCM_count; // normalize distribution such that the integral yields 1
	  double eCM = eCMmin +(j+0.5)*bin_width;
	  double gamma = 1.0+eCM/Mreduced; 
	  double vCM = Clight*sqrt(1.0-1.0/(gamma*gamma));  // with Clight from sigmarr.h
	  if (distri_flag)
	    { // output of energy distribution and binned DR cross section
	      fout << eCM << ", " << dist_eCM[j]/bin_width << ", " << sigmaDR[j] << endl;
	    }
	  else
	    { // integration v(eCM)*sigma(eCM)*f(eCM)*deCM, note that dist_eCM = f(eCM)*deCM
	      alphaDR += vCM*sigmaDR[j]*dist_eCM[j];
	      alphaRR += vCM/eCM*sigmarrscl(eCM,RRzeff,RRfraction,false,RRnmin,RRnmax,RRlmin)*dist_eCM[j];
	              // note that sigmarrscl comes with a factor eCM
	    }
	} // end for(j ...), end of convolution loop
      // *******************************************
      delete[] sigmaDR;
      delete[] dist_eCM;
      if (distri_flag) break;
      
      fout << eLAB << ",  " << erel << ", " << alphaDR << ", " << alphaRR << ", " << alphaDR+alphaRR << endl;
      if (!diagnostic_flag) { cout << "." << flush; }
    } // end for(n ...), end of loop over LAB energies
  // ***************************** cleanup ******************************
  fout.close();
  delete[] wDR;
  delete[] sDR;
  delete[] eDR;
  cout << endl << endl <<  " Data written to output file: " << outfilename << endl << endl;
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief measured angle in ESR cooler due to imperfections of magnetic guiding field
 *
 * @param zpos position in cooler along z-axis in cm (0 < zpos < 250)
 * @param angle on exit: angle at position zpos
 *
 * @return true, if zpos out of range
 * @return false, if zpos is ok
 */
bool CoolerAngleESR(double zpos, double &angle)
{
  // the angles have been corrected for an overall offset of  0.48 mrad
  const double theta[2501] =
{0.06216, 0.06145, 0.06075, 0.06004, 0.05934,
0.05863, 0.05792, 0.05722, 0.05651, 0.05581,
0.0551, 0.05448, 0.05385, 0.05323, 0.0526,
0.05198, 0.05135, 0.05072, 0.0501, 0.04947,
0.04885, 0.04829, 0.04774, 0.04718, 0.04663,
0.04607, 0.04552, 0.04496, 0.04441, 0.04385,
0.0433, 0.04281, 0.04232, 0.04183, 0.04134,
0.04084, 0.04035, 0.03986, 0.03937, 0.03888,
0.03839, 0.03795, 0.03752, 0.03708, 0.03665,
0.03621, 0.03577, 0.03534, 0.0349, 0.03447,
0.03403, 0.03364, 0.03326, 0.03287, 0.03249,
0.0321, 0.03171, 0.03133, 0.03094, 0.03056,
0.03017, 0.02983, 0.02948, 0.02914, 0.0288,
0.02846, 0.02811, 0.02777, 0.02743, 0.02708,
0.02674, 0.02644, 0.02613, 0.02583, 0.02553,
0.02522, 0.02492, 0.02462, 0.02432, 0.02401,
0.02371, 0.02344, 0.02317, 0.0229, 0.02263,
0.02236, 0.0221, 0.02183, 0.02156, 0.02129,
0.02102, 0.02078, 0.02054, 0.0203, 0.02006,
0.01983, 0.01959, 0.01935, 0.01911, 0.01887,
0.01863, 0.01842, 0.01821, 0.018, 0.01779,
0.01758, 0.01736, 0.01715, 0.01694, 0.01673,
0.01652, 0.01633, 0.01614, 0.01596, 0.01577,
0.01558, 0.01539, 0.0152, 0.01502, 0.01483,
0.01464, 0.01447, 0.01431, 0.01414, 0.01398,
0.01381, 0.01364, 0.01348, 0.01331, 0.01315,
0.01298, 0.01283, 0.01269, 0.01254, 0.01239,
0.01224, 0.0121, 0.01195, 0.0118, 0.01166,
0.01151, 0.01138, 0.01125, 0.01112, 0.01099,
0.01086, 0.01072, 0.01059, 0.01046, 0.01033,
0.0102, 0.01008, 0.00997, 0.00985, 0.00974,
0.00962, 0.0095, 0.00939, 0.00927, 0.00916,
0.00904, 0.00894, 0.00884, 0.00873, 0.00863,
0.00853, 0.00843, 0.00833, 0.00822, 0.00812,
0.00802, 0.00793, 0.00784, 0.00775, 0.00766,
0.00756, 0.00747, 0.00738, 0.00729, 0.0072,
0.00711, 0.00703, 0.00695, 0.00687, 0.00679,
0.00671, 0.00662, 0.00654, 0.00646, 0.00638,
0.0063, 0.00623, 0.00616, 0.00608, 0.00601,
0.00594, 0.00587, 0.0058, 0.00572, 0.00565,
0.00558, 0.00552, 0.00545, 0.00539, 0.00533,
0.00527, 0.0052, 0.00514, 0.00508, 0.00501,
0.00495, 0.00489, 0.00484, 0.00478, 0.00473,
0.00467, 0.00461, 0.00456, 0.0045, 0.00445,
0.00439, 0.00434, 0.00429, 0.00424, 0.00419,
0.00414, 0.00409, 0.00404, 0.00399, 0.00394,
0.00389, 0.00385, 0.0038, 0.00376, 0.00371,
0.00367, 0.00363, 0.00358, 0.00354, 0.00349,
0.00345, 0.00342, 0.00338, 0.00335, 0.00332,
0.00328, 0.00325, 0.00322, 0.00319, 0.00315,
0.00312, 0.00308, 0.00305, 0.00301, 0.00298,
0.00294, 0.0029, 0.00287, 0.00283, 0.0028,
0.00276, 0.00273, 0.0027, 0.00266, 0.00263,
0.0026, 0.00257, 0.00254, 0.0025, 0.00247,
0.00244, 0.00241, 0.00238, 0.00236, 0.00233,
0.0023, 0.00227, 0.00224, 0.00222, 0.00219,
0.00216, 0.00214, 0.00211, 0.00209, 0.00206,
0.00203, 0.00201, 0.00199, 0.00196, 0.00194,
0.00191, 0.00189, 0.00187, 0.00184, 0.00182,
0.0018, 0.00178, 0.00176, 0.00173, 0.00171,
0.00169, 0.00167, 0.00165, 0.00163, 0.00161,
0.0016, 0.00158, 0.00156, 0.00154, 0.00152,
0.0015, 0.00148, 0.00146, 0.00145, 0.00143,
0.00141, 0.00139, 0.00137, 0.00136, 0.00134,
0.00132, 0.00131, 0.00129, 0.00128, 0.00126,
0.00125, 0.00123, 0.00121, 0.0012, 0.00119,
0.00117, 0.00116, 0.00114, 0.00113, 0.00111,
0.0011, 0.00109, 0.00107, 0.00106, 0.00104,
0.00103, 0.00102, 0.00101, 0.000995741, 0.000984321,
0.000972901, 0.000961481, 0.000950061, 0.000938642, 0.000927222,
0.000915802, 0.000905462, 0.000895122, 0.000884781, 0.000874441,
0.000864101, 0.000853761, 0.000843421, 0.00083308, 0.00082274,
0.0008124, 0.000803041, 0.000793682, 0.000784323, 0.000774964,
0.000765606, 0.000756247, 0.000746888, 0.000737529, 0.00072817,
0.000718811, 0.000711031, 0.00070325, 0.00069547, 0.000687689,
0.000679909, 0.000672129, 0.000664348, 0.000656568, 0.000648787,
0.000641007, 0.000634271, 0.000627534, 0.000620798, 0.000614062,
0.000607326, 0.000600589, 0.000593853, 0.000587117, 0.00058038,
0.000573644, 0.000567285, 0.000560926, 0.000554567, 0.000548208,
0.000541849, 0.00053549, 0.000529131, 0.000522772, 0.000516413,
0.000510054, 0.000504092, 0.00049813, 0.000492168, 0.000486206,
0.000480244, 0.000474282, 0.00046832, 0.000462358, 0.000456396,
0.000450434, 0.000446974, 0.000443514, 0.000440053, 0.000436593,
0.000433133, 0.000429673, 0.000426213, 0.000422752, 0.000419292,
0.000415832, 0.000411735, 0.000407639, 0.000403542, 0.000399445,
0.000395348, 0.000391252, 0.000387155, 0.000383058, 0.000378962,
0.000374865, 0.000371337, 0.000367808, 0.00036428, 0.000360751,
0.000357223, 0.000353695, 0.000350166, 0.000346638, 0.000343109,
0.000339581, 0.000337572, 0.000335563, 0.000333553, 0.000331544,
0.000329535, 0.000327526, 0.000325517, 0.000323507, 0.000321498,
0.000319489, 0.000318455, 0.000317421, 0.000316387, 0.000315353,
0.000314319, 0.000313284, 0.00031225, 0.000311216, 0.000310182,
0.000309148, 0.000306429, 0.00030371, 0.000300991, 0.000298272,
0.000295552, 0.000292833, 0.000290114, 0.000287395, 0.000284676,
0.000281957, 0.000282048, 0.000282139, 0.00028223, 0.000282321,
0.000282413, 0.000282504, 0.000282595, 0.000282686, 0.000282777,
0.000282868, 0.000282918, 0.000282967, 0.000283016, 0.000283066,
0.000283115, 0.000283165, 0.000283215, 0.000283264, 0.000283314,
0.000283363, 0.000283326, 0.000283289, 0.000283253, 0.000283216,
0.000283179, 0.000283142, 0.000283105, 0.000283069, 0.000283032,
0.000282995, 0.000282928, 0.000282861, 0.000282795, 0.000282728,
0.000282661, 0.000282594, 0.000282527, 0.000282461, 0.000282394,
0.000282327, 0.000283347, 0.000284366, 0.000285386, 0.000286405,
0.000287425, 0.000288444, 0.000289464, 0.000290483, 0.000291503,
0.000292522, 0.000291682, 0.000290841, 0.000290001, 0.00028916,
0.00028832, 0.00028748, 0.000286639, 0.000285799, 0.000284958,
0.000284118, 0.000284701, 0.000285283, 0.000285866, 0.000286449,
0.000287032, 0.000287614, 0.000288197, 0.00028878, 0.000289362,
0.000289945, 0.000290553, 0.000291161, 0.000291769, 0.000292377,
0.000292984, 0.000293592, 0.0002942, 0.000294808, 0.000295416,
0.000296024, 0.000296657, 0.000297289, 0.000297922, 0.000298555,
0.000299187, 0.00029982, 0.000300453, 0.000301086, 0.000301718,
0.000302351, 0.000303532, 0.000304713, 0.000305894, 0.000307075,
0.000308257, 0.000309438, 0.000310619, 0.0003118, 0.000312981,
0.000314162, 0.000315679, 0.000317195, 0.000318712, 0.000320228,
0.000321745, 0.000323261, 0.000324778, 0.000326294, 0.000327811,
0.000329327, 0.000330194, 0.000331062, 0.000331929, 0.000332797,
0.000333664, 0.000334531, 0.000335399, 0.000336266, 0.000337134,
0.000338001, 0.000339075, 0.00034015, 0.000341224, 0.000342299,
0.000343373, 0.000344447, 0.000345522, 0.000346596, 0.000347671,
0.000348745, 0.000350017, 0.000351289, 0.000352561, 0.000353833,
0.000355105, 0.000356377, 0.000357649, 0.000358921, 0.000360193,
0.000361465, 0.000362366, 0.000363268, 0.000364169, 0.00036507,
0.000365972, 0.000366873, 0.000367774, 0.000368675, 0.000369577,
0.000370478, 0.000371106, 0.000371733, 0.000372361, 0.000372989,
0.000373617, 0.000374244, 0.000374872, 0.0003755, 0.000376127,
0.000376755, 0.00037699, 0.000377225, 0.00037746, 0.000377695,
0.000377931, 0.000378166, 0.000378401, 0.000378636, 0.000378871,
0.000379106, 0.000379303, 0.0003795, 0.000379698, 0.000379895,
0.000380092, 0.000380289, 0.000380486, 0.000380684, 0.000380881,
0.000381078, 0.000381259, 0.00038144, 0.00038162, 0.000381801,
0.000381982, 0.000382163, 0.000382344, 0.000382524, 0.000382705,
0.000382886, 0.000382704, 0.000382522, 0.00038234, 0.000382158,
0.000381977, 0.000381795, 0.000381613, 0.000381431, 0.000381249,
0.000381067, 0.000381154, 0.000381242, 0.000381329, 0.000381416,
0.000381504, 0.000381591, 0.000381678, 0.000381765, 0.000381853,
0.00038194, 0.000382784, 0.000383628, 0.000384473, 0.000385317,
0.000386161, 0.000387005, 0.000387849, 0.000388694, 0.000389538,
0.000390382, 0.000390399, 0.000390415, 0.000390432, 0.000390448,
0.000390465, 0.000390481, 0.000390498, 0.000390514, 0.000390531,
0.000390547, 0.000390991, 0.000391435, 0.000391879, 0.000392323,
0.000392767, 0.000393212, 0.000393656, 0.0003941, 0.000394544,
0.000394988, 0.000396004, 0.00039702, 0.000398035, 0.000399051,
0.000400067, 0.000401083, 0.000402099, 0.000403114, 0.00040413,
0.000405146, 0.000406156, 0.000407166, 0.000408177, 0.000409187,
0.000410197, 0.000411207, 0.000412217, 0.000413228, 0.000414238,
0.000415248, 0.000415928, 0.000416608, 0.000417289, 0.000417969,
0.000418649, 0.000419329, 0.000420009, 0.00042069, 0.00042137,
0.00042205, 0.000423291, 0.000424532, 0.000425773, 0.000427014,
0.000428256, 0.000429497, 0.000430738, 0.000431979, 0.00043322,
0.000434461, 0.000434986, 0.00043551, 0.000436035, 0.00043656,
0.000437085, 0.000437609, 0.000438134, 0.000438659, 0.000439183,
0.000439708, 0.000439665, 0.000439623, 0.00043958, 0.000439538,
0.000439495, 0.000439452, 0.00043941, 0.000439367, 0.000439325,
0.000439282, 0.00043967, 0.000440059, 0.000440447, 0.000440836,
0.000441224, 0.000441612, 0.000442001, 0.000442389, 0.000442778,
0.000443166, 0.000443157, 0.000443148, 0.000443138, 0.000443129,
0.00044312, 0.000443111, 0.000443102, 0.000443092, 0.000443083,
0.000443074, 0.000443173, 0.000443272, 0.000443371, 0.00044347,
0.000443569, 0.000443668, 0.000443767, 0.000443866, 0.000443965,
0.000444064, 0.000444131, 0.000444198, 0.000444264, 0.000444331,
0.000444398, 0.000444465, 0.000444532, 0.000444598, 0.000444665,
0.000444732, 0.000444783, 0.000444835, 0.000444886, 0.000444938,
0.000444989, 0.00044504, 0.000445092, 0.000445143, 0.000445195,
0.000445246, 0.000445545, 0.000445845, 0.000446144, 0.000446443,
0.000446743, 0.000447042, 0.000447341, 0.00044764, 0.00044794,
0.000448239, 0.000448746, 0.000449253, 0.000449761, 0.000450268,
0.000450775, 0.000451282, 0.000451789, 0.000452297, 0.000452804,
0.000453311, 0.000453799, 0.000454287, 0.000454775, 0.000455263,
0.000455751, 0.000456239, 0.000456727, 0.000457215, 0.000457703,
0.000458191, 0.000458927, 0.000459664, 0.000460401, 0.000461137,
0.000461874, 0.00046261, 0.000463347, 0.000464083, 0.00046482,
0.000465556, 0.000466107, 0.000466659, 0.00046721, 0.000467762,
0.000468313, 0.000468864, 0.000469416, 0.000469967, 0.000470519,
0.00047107, 0.000471506, 0.000471942, 0.000472378, 0.000472814,
0.00047325, 0.000473685, 0.000474121, 0.000474557, 0.000474993,
0.000475429, 0.000475697, 0.000475965, 0.000476232, 0.0004765,
0.000476768, 0.000477036, 0.000477304, 0.000477571, 0.000477839,
0.000478107, 0.000478097, 0.000478087, 0.000478078, 0.000478068,
0.000478058, 0.000478048, 0.000478038, 0.000478029, 0.000478019,
0.000478009, 0.000478112, 0.000478215, 0.000478317, 0.00047842,
0.000478523, 0.000478626, 0.000478729, 0.000478831, 0.000478934,
0.000479037, 0.000479305, 0.000479572, 0.00047984, 0.000480108,
0.000480376, 0.000480643, 0.000480911, 0.000481179, 0.000481446,
0.000481714, 0.000481755, 0.000481795, 0.000481836, 0.000481876,
0.000481917, 0.000481957, 0.000481998, 0.000482038, 0.000482079,
0.000482119, 0.000482206, 0.000482293, 0.000482379, 0.000482466,
0.000482553, 0.00048264, 0.000482727, 0.000482813, 0.0004829,
0.000482987, 0.000483052, 0.000483117, 0.000483182, 0.000483247,
0.000483312, 0.000483376, 0.000483441, 0.000483506, 0.000483571,
0.000483636, 0.000483845, 0.000484054, 0.000484263, 0.000484472,
0.000484681, 0.00048489, 0.000485099, 0.000485308, 0.000485517,
0.000485726, 0.00048624, 0.000486755, 0.000487269, 0.000487784,
0.000488298, 0.000488812, 0.000489327, 0.000489841, 0.000490356,
0.00049087, 0.00049129, 0.000491711, 0.000492131, 0.000492552,
0.000492972, 0.000493392, 0.000493813, 0.000494233, 0.000494654,
0.000495074, 0.000495396, 0.000495718, 0.000496041, 0.000496363,
0.000496685, 0.000497007, 0.000497329, 0.000497652, 0.000497974,
0.000498296, 0.000498734, 0.000499171, 0.000499609, 0.000500046,
0.000500484, 0.000500921, 0.000501358, 0.000501796, 0.000502234,
0.000502671, 0.000502891, 0.00050311, 0.00050333, 0.000503549,
0.000503768, 0.000503988, 0.000504207, 0.000504427, 0.000504646,
0.000504866, 0.000505106, 0.000505346, 0.000505585, 0.000505825,
0.000506065, 0.000506305, 0.000506545, 0.000506784, 0.000507024,
0.000507264, 0.000507319, 0.000507374, 0.000507429, 0.000507484,
0.000507539, 0.000507594, 0.000507649, 0.000507704, 0.000507759,
0.000507814, 0.000508189, 0.000508564, 0.00050894, 0.000509315,
0.00050969, 0.000510065, 0.00051044, 0.000510816, 0.000511191,
0.000511566, 0.000511841, 0.000512117, 0.000512392, 0.000512667,
0.000512943, 0.000513218, 0.000513493, 0.000513768, 0.000514044,
0.000514319, 0.000515365, 0.000516411, 0.000517457, 0.000518503,
0.000519549, 0.000520595, 0.000521641, 0.000522687, 0.000523733,
0.000524779, 0.000524736, 0.000524694, 0.000524651, 0.000524609,
0.000524566, 0.000524523, 0.000524481, 0.000524438, 0.000524396,
0.000524353, 0.000524421, 0.000524488, 0.000524556, 0.000524623,
0.000524691, 0.000524759, 0.000524826, 0.000524894, 0.000524961,
0.000525029, 0.000525076, 0.000525123, 0.000525171, 0.000525218,
0.000525265, 0.000525312, 0.000525359, 0.000525407, 0.000525454,
0.000525501, 0.000525978, 0.000526455, 0.000526931, 0.000527408,
0.000527885, 0.000528362, 0.000528839, 0.000529315, 0.000529792,
0.000530269, 0.000529511, 0.000528752, 0.000527994, 0.000527236,
0.000526477, 0.000525719, 0.000524961, 0.000524203, 0.000523444,
0.000522686, 0.000522719, 0.000522753, 0.000522786, 0.000522819,
0.000522853, 0.000522886, 0.000522919, 0.000522952, 0.000522986,
0.000523019, 0.000523272, 0.000523525, 0.000523778, 0.000524031,
0.000524285, 0.000524538, 0.000524791, 0.000525044, 0.000525297,
0.00052555, 0.000525733, 0.000525917, 0.0005261, 0.000526283,
0.000526467, 0.00052665, 0.000526833, 0.000527016, 0.0005272,
0.000527383, 0.000527496, 0.000527608, 0.000527721, 0.000527833,
0.000527946, 0.000528059, 0.000528171, 0.000528284, 0.000528396,
0.000528509, 0.000528299, 0.000528089, 0.000527879, 0.000527669,
0.000527459, 0.000527249, 0.000527039, 0.000526829, 0.000526619,
0.000526409, 0.000526598, 0.000526788, 0.000526977, 0.000527166,
0.000527356, 0.000527545, 0.000527734, 0.000527923, 0.000528113,
0.000528302, 0.00052832, 0.000528338, 0.000528357, 0.000528375,
0.000528393, 0.000528411, 0.000528429, 0.000528448, 0.000528466,
0.000528484, 0.000528522, 0.00052856, 0.000528598, 0.000528636,
0.000528674, 0.000528713, 0.000528751, 0.000528789, 0.000528827,
0.000528865, 0.000528603, 0.00052834, 0.000528078, 0.000527815,
0.000527553, 0.000527291, 0.000527028, 0.000526766, 0.000526503,
0.000526241, 0.00052686, 0.000527478, 0.000528097, 0.000528715,
0.000529334, 0.000529953, 0.000530571, 0.00053119, 0.000531808,
0.000532427, 0.000532679, 0.000532931, 0.000533182, 0.000533434,
0.000533686, 0.000533938, 0.00053419, 0.000534441, 0.000534693,
0.000534945, 0.000535001, 0.000535056, 0.000535111, 0.000535167,
0.000535223, 0.000535278, 0.000535333, 0.000535389, 0.000535445,
0.0005355, 0.000535288, 0.000535076, 0.000534864, 0.000534652,
0.00053444, 0.000534228, 0.000534016, 0.000533804, 0.000533592,
0.00053338, 0.000533182, 0.000532983, 0.000532785, 0.000532586,
0.000532388, 0.000532189, 0.000531991, 0.000531792, 0.000531594,
0.000531395, 0.000531162, 0.000530928, 0.000530695, 0.000530462,
0.000530228, 0.000529995, 0.000529762, 0.000529529, 0.000529295,
0.000529062, 0.000528601, 0.000528139, 0.000527678, 0.000527217,
0.000526756, 0.000526294, 0.000525833, 0.000525372, 0.00052491,
0.000524449, 0.000523872, 0.000523296, 0.000522719, 0.000522142,
0.000521566, 0.000520989, 0.000520412, 0.000519835, 0.000519259,
0.000518682, 0.00051864, 0.000518599, 0.000518557, 0.000518516,
0.000518474, 0.000518432, 0.000518391, 0.000518349, 0.000518308,
0.000518266, 0.000517943, 0.000517619, 0.000517296, 0.000516973,
0.00051665, 0.000516326, 0.000516003, 0.00051568, 0.000515356,
0.000515033, 0.000514652, 0.00051427, 0.000513889, 0.000513507,
0.000513126, 0.000512744, 0.000512362, 0.000511981, 0.0005116,
0.000511218, 0.000511405, 0.000511593, 0.00051178, 0.000511967,
0.000512155, 0.000512342, 0.000512529, 0.000512716, 0.000512904,
0.000513091, 0.000513156, 0.000513221, 0.000513286, 0.000513351,
0.000513416, 0.000513481, 0.000513546, 0.000513611, 0.000513676,
0.000513741, 0.000513685, 0.000513629, 0.000513573, 0.000513517,
0.000513461, 0.000513405, 0.000513349, 0.000513293, 0.000513237,
0.000513181, 0.000513419, 0.000513657, 0.000513894, 0.000514132,
0.00051437, 0.000514608, 0.000514846, 0.000515083, 0.000515321,
0.000515559, 0.000515797, 0.000516035, 0.000516274, 0.000516512,
0.00051675, 0.000516988, 0.000517226, 0.000517465, 0.000517703,
0.000517941, 0.000517974, 0.000518007, 0.00051804, 0.000518073,
0.000518107, 0.00051814, 0.000518173, 0.000518206, 0.000518239,
0.000518272, 0.000519111, 0.00051995, 0.000520789, 0.000521628,
0.000522467, 0.000523305, 0.000524144, 0.000524983, 0.000525822,
0.000526661, 0.000527174, 0.000527687, 0.0005282, 0.000528713,
0.000529225, 0.000529738, 0.000530251, 0.000530764, 0.000531277,
0.00053179, 0.000532255, 0.00053272, 0.000533186, 0.000533651,
0.000534116, 0.000534581, 0.000535046, 0.000535512, 0.000535977,
0.000536442, 0.000536616, 0.00053679, 0.000536964, 0.000537138,
0.000537313, 0.000537487, 0.000537661, 0.000537835, 0.000538009,
0.000538183, 0.000538034, 0.000537885, 0.000537736, 0.000537587,
0.000537438, 0.000537288, 0.000537139, 0.00053699, 0.000536841,
0.000536692, 0.000536108, 0.000535524, 0.00053494, 0.000534356,
0.000533772, 0.000533187, 0.000532603, 0.000532019, 0.000531435,
0.000530851, 0.000530599, 0.000530347, 0.000530094, 0.000529842,
0.00052959, 0.000529338, 0.000529086, 0.000528833, 0.000528581,
0.000528329, 0.000527972, 0.000527616, 0.00052726, 0.000526903,
0.000526546, 0.00052619, 0.000525833, 0.000525477, 0.00052512,
0.000524764, 0.00052466, 0.000524557, 0.000524453, 0.000524349,
0.000524245, 0.000524142, 0.000524038, 0.000523934, 0.000523831,
0.000523727, 0.000523255, 0.000522782, 0.00052231, 0.000521838,
0.000521365, 0.000520893, 0.000520421, 0.000519949, 0.000519476,
0.000519004, 0.000518694, 0.000518383, 0.000518073, 0.000517762,
0.000517452, 0.000517142, 0.000516831, 0.000516521, 0.00051621,
0.0005159, 0.000516428, 0.000516956, 0.000517484, 0.000518012,
0.00051854, 0.000519067, 0.000519595, 0.000520123, 0.000520651,
0.000521179, 0.000521662, 0.000522145, 0.000522628, 0.000523111,
0.000523593, 0.000524076, 0.000524559, 0.000525042, 0.000525525,
0.000526008, 0.000526765, 0.000527523, 0.00052828, 0.000529038,
0.000529795, 0.000530552, 0.00053131, 0.000532067, 0.000532825,
0.000533582, 0.000534644, 0.000535707, 0.000536769, 0.000537831,
0.000538894, 0.000539956, 0.000541018, 0.00054208, 0.000543143,
0.000544205, 0.00054402, 0.000543835, 0.00054365, 0.000543465,
0.00054328, 0.000543095, 0.00054291, 0.000542725, 0.00054254,
0.000542355, 0.000542067, 0.00054178, 0.000541493, 0.000541205,
0.000540918, 0.00054063, 0.000540343, 0.000540055, 0.000539768,
0.00053948, 0.000540022, 0.000540564, 0.000541106, 0.000541648,
0.000542191, 0.000542733, 0.000543275, 0.000543817, 0.000544359,
0.000544901, 0.00054451, 0.000544119, 0.000543729, 0.000543338,
0.000542947, 0.000542556, 0.000542165, 0.000541775, 0.000541384,
0.000540993, 0.000540454, 0.000539915, 0.000539375, 0.000538836,
0.000538297, 0.000537758, 0.000537219, 0.000536679, 0.00053614,
0.000535601, 0.000535888, 0.000536175, 0.000536463, 0.00053675,
0.000537037, 0.000537324, 0.000537611, 0.000537899, 0.000538186,
0.000538473, 0.000538137, 0.0005378, 0.000537464, 0.000537128,
0.000536791, 0.000536455, 0.000536119, 0.000535783, 0.000535446,
0.00053511, 0.000534523, 0.000533935, 0.000533348, 0.000532761,
0.000532174, 0.000531586, 0.000530999, 0.000530412, 0.000529824,
0.000529237, 0.000529089, 0.000528942, 0.000528794, 0.000528646,
0.000528499, 0.000528351, 0.000528203, 0.000528055, 0.000527908,
0.00052776, 0.000528237, 0.000528713, 0.00052919, 0.000529666,
0.000530143, 0.00053062, 0.000531096, 0.000531573, 0.000532049,
0.000532526, 0.000533038, 0.00053355, 0.000534062, 0.000534574,
0.000535086, 0.000535598, 0.00053611, 0.000536622, 0.000537134,
0.000537646, 0.000538005, 0.000538363, 0.000538722, 0.00053908,
0.000539439, 0.000539798, 0.000540156, 0.000540515, 0.000540873,
0.000541232, 0.000541044, 0.000540856, 0.000540667, 0.000540479,
0.000540291, 0.000540103, 0.000539915, 0.000539726, 0.000539538,
0.00053935, 0.000538807, 0.000538264, 0.00053772, 0.000537177,
0.000536634, 0.000536091, 0.000535548, 0.000535004, 0.000534461,
0.000533918, 0.000533905, 0.000533893, 0.00053388, 0.000533867,
0.000533855, 0.000533842, 0.000533829, 0.000533816, 0.000533804,
0.000533791, 0.000533486, 0.000533181, 0.000532877, 0.000532572,
0.000532267, 0.000531962, 0.000531657, 0.000531353, 0.000531048,
0.000530743, 0.000530166, 0.00052959, 0.000529013, 0.000528437,
0.00052786, 0.000527283, 0.000526707, 0.00052613, 0.000525554,
0.000524977, 0.000524357, 0.000523738, 0.000523118, 0.000522498,
0.000521879, 0.000521259, 0.000520639, 0.000520019, 0.0005194,
0.00051878, 0.000518603, 0.000518426, 0.000518249, 0.000518072,
0.000517895, 0.000517718, 0.000517541, 0.000517364, 0.000517187,
0.00051701, 0.000516063, 0.000515116, 0.000514169, 0.000513222,
0.000512276, 0.000511329, 0.000510382, 0.000509435, 0.000508488,
0.000507541, 0.000507323, 0.000507105, 0.000506886, 0.000506668,
0.00050645, 0.000506232, 0.000506014, 0.000505795, 0.000505577,
0.000505359, 0.000505332, 0.000505305, 0.000505278, 0.000505251,
0.000505225, 0.000505198, 0.000505171, 0.000505144, 0.000505117,
0.00050509, 0.000505341, 0.000505593, 0.000505844, 0.000506096,
0.000506347, 0.000506598, 0.00050685, 0.000507101, 0.000507353,
0.000507604, 0.000508071, 0.000508539, 0.000509006, 0.000509473,
0.000509941, 0.000510408, 0.000510875, 0.000511342, 0.00051181,
0.000512277, 0.000512978, 0.000513679, 0.000514379, 0.00051508,
0.000515781, 0.000516482, 0.000517183, 0.000517883, 0.000518584,
0.000519285, 0.000519377, 0.000519469, 0.000519561, 0.000519653,
0.000519745, 0.000519837, 0.000519929, 0.000520021, 0.000520113,
0.000520205, 0.000520173, 0.000520142, 0.00052011, 0.000520078,
0.000520047, 0.000520015, 0.000519983, 0.000519951, 0.00051992,
0.000519888, 0.000519905, 0.000519922, 0.000519939, 0.000519956,
0.000519973, 0.00051999, 0.000520007, 0.000520024, 0.000520041,
0.000520058, 0.000519915, 0.000519772, 0.000519629, 0.000519486,
0.000519343, 0.000519201, 0.000519058, 0.000518915, 0.000518772,
0.000518629, 0.000518054, 0.000517479, 0.000516903, 0.000516328,
0.000515753, 0.000515178, 0.000514603, 0.000514027, 0.000513452,
0.000512877, 0.000512773, 0.00051267, 0.000512566, 0.000512462,
0.000512359, 0.000512255, 0.000512151, 0.000512047, 0.000511944,
0.00051184, 0.000511824, 0.000511807, 0.000511791, 0.000511775,
0.000511759, 0.000511742, 0.000511726, 0.00051171, 0.000511693,
0.000511677, 0.000512367, 0.000513057, 0.000513747, 0.000514437,
0.000515127, 0.000515818, 0.000516508, 0.000517198, 0.000517888,
0.000518578, 0.000518656, 0.000518734, 0.000518812, 0.00051889,
0.000518968, 0.000519046, 0.000519124, 0.000519202, 0.00051928,
0.000519358, 0.000520128, 0.000520899, 0.000521669, 0.00052244,
0.00052321, 0.00052398, 0.000524751, 0.000525521, 0.000526292,
0.000527062, 0.000527607, 0.000528153, 0.000528698, 0.000529244,
0.000529789, 0.000530334, 0.00053088, 0.000531425, 0.000531971,
0.000532516, 0.000533411, 0.000534305, 0.0005352, 0.000536094,
0.000536989, 0.000537883, 0.000538778, 0.000539672, 0.000540566,
0.000541461, 0.000541601, 0.000541741, 0.000541881, 0.000542021,
0.000542162, 0.000542302, 0.000542442, 0.000542582, 0.000542722,
0.000542862, 0.000544118, 0.000545375, 0.000546631, 0.000547887,
0.000549144, 0.0005504, 0.000551656, 0.000552912, 0.000554169,
0.000555425, 0.000556221, 0.000557016, 0.000557812, 0.000558608,
0.000559404, 0.000560199, 0.000560995, 0.000561791, 0.000562586,
0.000563382, 0.000564182, 0.000564983, 0.000565783, 0.000566583,
0.000567384, 0.000568184, 0.000568984, 0.000569784, 0.000570585,
0.000571385, 0.000572218, 0.000573051, 0.000573884, 0.000574717,
0.00057555, 0.000576384, 0.000577217, 0.00057805, 0.000578883,
0.000579716, 0.00058035, 0.000580985, 0.000581619, 0.000582254,
0.000582889, 0.000583523, 0.000584158, 0.000584792, 0.000585427,
0.000586061, 0.000586728, 0.000587394, 0.000588061, 0.000588728,
0.000589394, 0.000590061, 0.000590728, 0.000591395, 0.000592061,
0.000592728, 0.000593476, 0.000594225, 0.000594973, 0.000595722,
0.00059647, 0.000597218, 0.000597967, 0.000598715, 0.000599464,
0.000600212, 0.000601466, 0.00060272, 0.000603973, 0.000605227,
0.000606481, 0.000607735, 0.000608989, 0.000610242, 0.000611496,
0.00061275, 0.000614047, 0.000615343, 0.00061664, 0.000617936,
0.000619233, 0.000620529, 0.000621826, 0.000623122, 0.000624419,
0.000625715, 0.000627924, 0.000630132, 0.000632341, 0.00063455,
0.000636759, 0.000638967, 0.000641176, 0.000643385, 0.000645593,
0.000647802, 0.00064976, 0.000651718, 0.000653676, 0.000655634,
0.000657592, 0.000659549, 0.000661507, 0.000663465, 0.000665423,
0.000667381, 0.000670074, 0.000672767, 0.00067546, 0.000678153,
0.000680846, 0.000683539, 0.000686232, 0.000688925, 0.000691618,
0.000694311, 0.000697682, 0.000701052, 0.000704423, 0.000707794,
0.000711164, 0.000714535, 0.000717906, 0.000721277, 0.000724647,
0.000728018, 0.000731152, 0.000734287, 0.000737421, 0.000740556,
0.000743691, 0.000746825, 0.00074996, 0.000753094, 0.000756229,
0.000759363, 0.000762827, 0.000766291, 0.000769755, 0.000773219,
0.000776683, 0.000780147, 0.000783611, 0.000787075, 0.000790539,
0.000794003, 0.000798409, 0.000802815, 0.000807222, 0.000811628,
0.000816034, 0.00082044, 0.000824846, 0.000829253, 0.000833659,
0.000838065, 0.000842554, 0.000847044, 0.000851533, 0.000856022,
0.000860511, 0.000865001, 0.00086949, 0.000873979, 0.000878469,
0.000882958, 0.000887472, 0.000891987, 0.000896501, 0.000901016,
0.00090553, 0.000910044, 0.000914559, 0.000919073, 0.000923588,
0.000928102, 0.000933822, 0.000939542, 0.000945261, 0.000950981,
0.000956701, 0.000962421, 0.000968141, 0.00097386, 0.00097958,
0.0009853, 0.00099177, 0.00099824, 0.001, 0.00101,
0.00102, 0.00102, 0.00103, 0.00104, 0.00104,
0.00105, 0.00106, 0.00106, 0.00107, 0.00108,
0.00109, 0.00109, 0.0011, 0.00111, 0.00111,
0.00112, 0.00113, 0.00114, 0.00114, 0.00115,
0.00116, 0.00117, 0.00118, 0.00118, 0.00119,
0.0012, 0.00121, 0.00122, 0.00123, 0.00124,
0.00125, 0.00126, 0.00127, 0.00128, 0.00129,
0.0013, 0.00131, 0.00132, 0.00133, 0.00134,
0.00135, 0.00137, 0.00138, 0.00139, 0.0014,
0.00141, 0.00142, 0.00144, 0.00145, 0.00146,
0.00148, 0.00149, 0.0015, 0.00151, 0.00153,
0.00154, 0.00155, 0.00157, 0.00158, 0.0016,
0.00161, 0.00162, 0.00164, 0.00165, 0.00167,
0.00168, 0.0017, 0.00171, 0.00173, 0.00174,
0.00176, 0.00177, 0.00179, 0.0018, 0.00182,
0.00183, 0.00185, 0.00186, 0.00188, 0.0019,
0.00192, 0.00193, 0.00195, 0.00197, 0.00198,
0.002, 0.00202, 0.00204, 0.00206, 0.00208,
0.0021, 0.00212, 0.00214, 0.00216, 0.00218,
0.0022, 0.00222, 0.00224, 0.00227, 0.00229,
0.00231, 0.00233, 0.00235, 0.00238, 0.0024,
0.00242, 0.00244, 0.00247, 0.00249, 0.00252,
0.00255, 0.00257, 0.0026, 0.00262, 0.00265,
0.00267, 0.0027, 0.00273, 0.00276, 0.00279,
0.00281, 0.00284, 0.00287, 0.0029, 0.00293,
0.00296, 0.00299, 0.00302, 0.00306, 0.00309,
0.00312, 0.00315, 0.00318, 0.00322, 0.00325,
0.00328, 0.00332, 0.00335, 0.00339, 0.00343,
0.00347, 0.0035, 0.00354, 0.00358, 0.00361,
0.00365, 0.00369, 0.00373, 0.00377, 0.00381,
0.00385, 0.00389, 0.00393, 0.00397, 0.00401,
0.00405, 0.00411, 0.00416, 0.00422, 0.00427,
0.00433, 0.00439, 0.00444, 0.0045, 0.00455,
0.00461, 0.00466, 0.00471, 0.00477, 0.00482,
0.00487, 0.00492, 0.00497, 0.00503, 0.00508,
0.00513, 0.00519, 0.00524, 0.0053, 0.00536,
0.00541, 0.00547, 0.00553, 0.00559, 0.00564,
0.0057, 0.00576, 0.00583, 0.00589, 0.00596,
0.00602, 0.00608, 0.00615, 0.00621, 0.00628,
0.00634, 0.00641, 0.00648, 0.00655, 0.00662,
0.00669, 0.00677, 0.00684, 0.00691, 0.00698,
0.00705, 0.00713, 0.00721, 0.00729, 0.00737,
0.00745, 0.00752, 0.0076, 0.00768, 0.00776,
0.00784, 0.00793, 0.00802, 0.00811, 0.0082,
0.00829, 0.00837, 0.00846, 0.00855, 0.00864,
0.00873, 0.00883, 0.00892, 0.00902, 0.00912,
0.00922, 0.00931, 0.00941, 0.00951, 0.0096,
0.0097, 0.00981, 0.00992, 0.01003, 0.01014,
0.01025, 0.01035, 0.01046, 0.01057, 0.01068,
0.01079, 0.01091, 0.01103, 0.01116, 0.01128,
0.0114, 0.01152, 0.01164, 0.01177, 0.01189,
0.01201, 0.01214, 0.01228, 0.01241, 0.01255,
0.01268, 0.01281, 0.01295, 0.01308, 0.01322,
0.01335, 0.0135, 0.01365, 0.0138, 0.01395,
0.0141, 0.01425, 0.0144, 0.01455, 0.0147,
0.01485, 0.01502, 0.01518, 0.01535, 0.01552,
0.01569, 0.01585, 0.01602, 0.01619, 0.01635,
0.01652, 0.01671, 0.01689, 0.01708, 0.01726,
0.01745, 0.01763, 0.01782, 0.018, 0.01819,
0.01837, 0.01858, 0.01878, 0.01899, 0.0192,
0.0194, 0.01961, 0.01982, 0.02003, 0.02023,
0.02044, 0.02067, 0.0209, 0.02113, 0.02136,
0.02159, 0.02181, 0.02204, 0.02227, 0.0225,
0.02273, 0.02298, 0.02324, 0.0235, 0.02375,
0.024, 0.02426, 0.02451, 0.02477, 0.02502,
0.02528, 0.02556, 0.02585, 0.02613, 0.02642,
0.0267, 0.02698, 0.02727, 0.02755, 0.02784,
0.02812, 0.02844, 0.02875, 0.02907, 0.02938,
0.0297, 0.03002, 0.03033, 0.03065, 0.03096,
0.03128, 0.03163, 0.03198, 0.03233, 0.03268,
0.03304, 0.03339, 0.03374, 0.03409, 0.03444,
0.03479, 0.03518, 0.03557, 0.03596, 0.03635,
0.03674, 0.03713, 0.03752, 0.03791, 0.0383,
0.03869, 0.03912, 0.03956, 0.03999, 0.04043,
0.04086, 0.04129, 0.04173, 0.04216, 0.0426,
0.04303, 0.04351, 0.044, 0.04448, 0.04496,
0.04544, 0.04593, 0.04641, 0.04689, 0.04738,
0.04786, 0.0484, 0.04894, 0.04947, 0.05001,
0.05055, 0.05109, 0.05163, 0.05216, 0.0527,
0.05324, 0.05384, 0.05443, 0.05503, 0.05563,
0.05622, 0.05682, 0.05742, 0.05802, 0.05861,
0.05921};
  
  double zmm = 10.0*zpos; // position in mm
  int i = floor(zmm);
  if ((i<0) || (i>2500)) return false;  // out of range
  
  angle = theta[i]+(theta[i+1]-theta[i])*(zmm-(double)i); // interpolation
  return false; // normal exit
}
