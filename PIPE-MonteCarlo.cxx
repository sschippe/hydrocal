/**
 * @file PIPE-MonteCarlo.cxx
 *
 * @brief Monte-Carlo simulation of molecular break-up in the PIPE experimental setup
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: PIPE-MonteCarlo.cxx 650 2019-09-25 07:30:07Z iamp $
   @endverbatim
 *
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cstring>
#include <chrono>      // clocks and time
#include <random>      // random number generators and distributions  
#include <functional>  // for "bind (random_generator, random_distribution)"
#include <eigen3/Eigen/Dense> // requires package libeigen3-dev
#include "svnrevision.h"

using namespace Eigen;
using namespace std;


typedef Matrix<float, 5, 5>  Matrix5f;
typedef Matrix<float, 5, 1>  Vector5f;
typedef Matrix<float, 1, 5>  RowVector5f;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief random polar angle
 *
 * @param r random number from a uniform distribution of real numbers from the interval [0,1]
 * @param a2 ansisotropy parameter of the distribution (allowed values: -1.0, 0.0, 1.0, 2.0)
 *
 * @return random value for the polar angle 
 */ 
double random_theta(double r, double &a2)
{
  // note that in general the angle theta is calculated as a root of the cubic equation
  // r = 0.25*[2-(a2-2)*cos(theta)+a2*cos(theta)^3].
  // Here, analytically solutions are provided for some special cases of a2.
  // The formulas have been derived by using Mathematica.
  // see also Amitay et al. PRA 54 (1996) 4032
  const double eps = 1.0E-9;
  const double Pi = 3.1415926535;
  const complex<double> I(0.0,1.0);
  const complex<double> minusonetwothirds(-0.5,0.866025403784439);
  double theta=0.0;
  if (fabs(a2)<eps)
    { // isotropic distribution
      a2 = 0.0;
      theta = acos(1.0-2.0*r);
    }
  else if (fabs(a2+1.0)<eps)
    {
      a2 = -1.0;
      complex<double> sqrtr=sqrt((1.0-r)*r);
      double aterm=arg(1.0-2.0*r+2.0*I*sqrtr)/3.0;
      theta = acos(sqrt(3.0)*sin(aterm)-cos(aterm));
    }
  else if (fabs(a2-1.0)<eps)
    { // does not work properly 
      a2 = 1.0;
      complex<double> sqrtr = sqrt(7.0/27.0-(1.0-r)*r);
      complex<double> term1 = exp(log(2.0*r-1.0-2.0*sqrtr)/3.0);
      complex<double> term2 = exp(log(1.0-2.0*r-2.0*sqrtr)/3.0);
      complex<double> term3 = minusonetwothirds*(term1-term2);
      theta = acos(real(term3));
    }
  else if (fabs(a2-2.0)<eps)
    {
      a2 = 2.0;
      double term = (r<0.5) ? -pow(1.0-2.0*r, 1.0/3.0) : pow(2.0*r-1.0, 1.0/3.0);
      theta = acos(term);
    }
  else if (fabs(a2-3.0)<eps)
    { // for testing puposes only
      a2 = 3.0;
      theta = r*Pi;
    }
  else
    { // isotropic distribution
      a2 = 0.0;
      theta = acos(1.0-2.0*r);
    }
  return theta;
}

Matrix5f matrix_drift(double L)
{
  Matrix5f m;
  m << 1, L, 0, 0, 0, 
    0, 1, 0, 0, 0,
    0, 0 ,1, 0, 0,
    0, 0, 0, 1, L, 
    0, 0, 0, 0, 1;
  return m;
}

Matrix5f matrix_lens(double f)
{
  Matrix5f m;
  m << 1, 0, 0, 0, 0, 
    1/f, 1, 0, 0, 0,
    0, 0 ,1, 0, 0,
    0, 0, 0, 1, 0, 
    0, 0, 0, -1/f, 1;
  return m;
}

Matrix5f matrix_dipole_magnet(double rho0, double phi0)
{ // see Wollnik, Eqs. 4.8a, 4.8b (also 4.25 for focussing)
  double c = cos(phi0);
  double s = sin(phi0);
  double d = rho0*(1.0-c);
  Matrix5f m;
  m << c, s*rho0, d, 0, 0,
    -s/rho0, c, s, 0, 0, 
     0, 0, 1, 0, 0, 
    0, 0, 0, 1, rho0*phi0, 
    0, 0, 0, 0, 1;
  return m;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Entry to Monte Carlo simulation of molecular breakup in the PIPE setup
 */
int PIPE_MonteCarlo()
{
  const double Pi=3.1415926535;
  
  const double PIPE_interact_length = 490;     // length (in mm) of interaction region
  const double PIPE_drift_demerger  = 865;     // drift length (in mm) from exit of interaction region to entrance of demerger
  const double PIPE_drift_detA      = 150;     // drift length (in mm) from exit of interaction region to detector A
  const double PIPE_drift_detB      = 500;     // drift length (in mm) from exit of demerger to detector B
  const double PIPE_dia_interact    = 6;       // diameter (in mm) of the aperture at the exit of the interaction region
  const double PIPE_dia_demerger    = 27;      // diameter (in mm) of the aperture at the entrance of the demerger 
  const double PIPE_dia_detectorA   = 75;      // diameter (in mm) of detector A
  const double PIPE_demerger_radius = 608.28;  // bending radius (in mm) of demerger magnet

  
  cout<< endl << " Monte Carlo simulation of molecular break up in the photo-ion merged-beams setup PIPE";
  cout<< endl;
  double interact_length = PIPE_interact_length;
  double drift_demerger  = PIPE_drift_demerger;
  double drift_detA      = PIPE_drift_detA;
  double drift_detB      = PIPE_drift_detB;
  double dia_interact    = PIPE_dia_interact;
  double dia_demerger    = PIPE_dia_demerger;
  double dia_detectorA   = PIPE_dia_detectorA;
  char answer;

  cout << endl << " Change default geometry values? (y/n) .................: ";
  cin >> answer;
  if ((answer=='y')||(answer=='Y'))
    {
      cout << endl << " Give length of interaction region (IR) in mm ..........: ";
      cin >> interact_length;
      cout << endl << " Give drift length from IR to detector A in mm .........: ";
      cin >> drift_detA;
      cout << endl << " Give outer diameter of Detector A in mm ...............: ";
      cin >> dia_detectorA;
      cout << endl << " Give drift length from IR to demerger (DEM) in mm .....: ";
      cin >> drift_demerger;
      cout << endl << " Give dia of aperture at entrance of DEM ...............: ";
      cin >> dia_demerger;
      cout << endl << " Give drift length from DEM to detector B in mm ........: ";
      cin >> drift_detB;
      cout << endl << " Give diameter of diaphragm at end of IR in mm .........: ";
      cin >> dia_interact; 
    }
  
  double ion_mass, ion_beam_fwhm, ion_beam_energy, frag1_mass, frag2_mass, frag_KER, frag_anisotropy=0.0;
  cout << endl << " Give primary ion energy in keV ........................: ";
  cin >> ion_beam_energy;
  ion_beam_energy *= 1000.0; // conversion from keV to eV
  cout << endl << " Give FWHM of ion beam in mm ...........................: ";
  cin >> ion_beam_fwhm;
  cout << endl << " Give mass of fragment 1 in u ..........................: ";
  cin >> frag1_mass;
  cout << endl << " Give mass of fragment 2 in u ..........................: ";
  cin >> frag2_mass;
  cout << endl << " Give kinetic energy release in eV .....................: ";
  cin >> frag_KER;
  cout << endl << " Give anisotropy parameter of fragment distribution ....: ";
  cin >> frag_anisotropy;


  // **** kinematic factors ****
  // initial momentum of primary molecular beam: p0=sqrt[2*Eion*(m1+m2)]
  double initial_momentum = sqrt(2.0*ion_beam_energy*(frag1_mass+frag2_mass));
  double initial_velocity = initial_momentum/(frag1_mass+frag2_mass);
  // momentum conservation upon fragmentation in comoving frame: p1=-p2=|pf|
  // kinetic energy release: KER=p1*p1/(2*m1)+p2*p2/(2*m2) = pf*pf/2*(m1+m2)/(m1*m2) => |pf| = sqrt[2*KER*m1*m2/(m1+m2)]
  double frag_momentum = sqrt(2.0*frag_KER*frag1_mass*frag2_mass/(frag1_mass+frag2_mass));
  double frag_momentum_ratio = frag_momentum/initial_momentum;
  double frag1_velocity = frag_momentum/frag1_mass;		   
  double frag2_velocity = -frag_momentum/frag2_mass;		   
  // distance from position of fragmentation event to detector: L
  // distance of particle 1 from axis: d1=v1*t*sin(theta)=v1*(L/v0)*sin(theta)
  // distance of particle 2 from axis: d2=v2*t*sin(theta)=v1*(L/v0)*sin(theta)
  double frag1_velocity_ratio = frag1_velocity/initial_velocity;
  double frag2_velocity_ratio = frag2_velocity/initial_velocity;


  double dia_trinos_x, dia_trinos_y;
  cout << endl << " Give diameter of x-Trinos slit at end of IR in mm .....: ";
  cin >> dia_trinos_x;
  cout << endl << " Give diameter of y-Trinos slit at end of IR in mm .....: ";
  cin >> dia_trinos_y;

  double rho = PIPE_demerger_radius;
  double f_demerger = rho/tan(Pi*26.565/180); // focal length due to oblique entrance and exit with an angle of 26.565 deg (double focussing)
  Matrix5f m_demerger = matrix_lens(f_demerger)*matrix_dipole_magnet(rho,0.5*Pi)*matrix_lens(f_demerger);
  Matrix5f m_drift_demerger = matrix_drift(drift_demerger);
  Matrix5f m_drift_detA = matrix_drift(drift_detA);
  Matrix5f m_drift_detB = matrix_drift(drift_detB);
  
  
  // ******************  seed for the  random generator  ***************************
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
 
  int MonteCarlo_iterations;
  cout << endl << " Give the number of Monte-Carlo iterations .............: ";
  cin >> MonteCarlo_iterations;
  
  // ********************* open output file and write header ************************         
  fstream fout;
  char outfilename[200], extension[5];
  sprintf(extension,".pmc");
  cout << endl << " Give filename for output (" << extension << ") .......................: ";
  cin >> outfilename;
  strcat(outfilename,extension);
  
  fout.open(outfilename,fstream::out);
  
  fout << "#######################################################################################################################" << endl;
  fout << "###  Monte Carlo simulation of molecular breakup in the PIPE setup" << endl;
  fout << "###" << endl;
  fout << "###  hydrocal SVN revision : " << SVNrevision << endl; // defined in SVNrevision.h
  fout << "###               filename : " << outfilename << endl;
  fout << "###  interact. length (mm) : " << interact_length << endl;
  fout << "###  apert. interact. (mm) : " << dia_interact << endl;
  fout << "###  apert. demerger  (mm) : " << dia_demerger << endl;
  fout << "###     x Trinos slit (mm) : " << dia_trinos_x << endl;
  fout << "###     y Trinos slit (mm) : " << dia_trinos_y << endl;
  fout << "### drift length detA (mm) : " << drift_detA << endl;
  fout << "###  diameter of detA (mm) : " << dia_detectorA << endl;
  fout << "### drift length dem  (mm) : " << drift_demerger << endl;
  fout << "###      demerger rho (mm) : " << rho << endl;
  fout << "### drift length detB (mm) : " << drift_detB << endl;
  fout << "###  ion beam energy (keV) : " << ion_beam_energy*1E-3 << endl;
  fout << "###     ion beam FWHM (mm) : " << ion_beam_fwhm << endl;
  fout << "###    fragment 1 mass (u) : " << frag1_mass << endl;
  fout << "###      fragment 1 qB (G) : " << 45.529*sqrt(frag1_mass*ion_beam_energy*1E-3)/(rho*0.001) << endl;  
  fout << "###    fragment 2 mass (u) : " << frag2_mass << endl;
  fout << "###      fragment 2 qB (G) : " << 45.529*sqrt(frag2_mass*ion_beam_energy*1E-3)/(rho*0.001) << endl;
  fout << "###               KER (eV) : " << frag_KER << endl;
  fout << "###   anisotropy parameter : " << frag_anisotropy << endl;
  fout << "###            random seed : " << seed << endl;
  fout << "### Monte-Carlo iterations : " << MonteCarlo_iterations << endl;
  fout << "#############################################################################################################################" << endl;
  fout << "###  detAxfrag1 (mm) [1]   detAyfrag1 (mm) [2]   detBxfrag1 (mm) [3]   detBdeltafrag1 [4]";
  fout << "   detAxfrag2 (mm) [5]   detAyfrag2 (mm) [6]   detBxfrag2 (mm) [7]   detBdeltafrag2 [8]" << endl;
  fout << "###--------------------------------------------------------------------------------------------------------------------------" << endl;
  
  cout << endl << endl << " now calculating, be patient " << flush;

  //************************ initialize random generators and random distributions  *******************

  mt19937_64 randgen(seed);  // 64bit "Mersenne Twister 19937" random number generator
  const double ln2=log(2.0);

  uniform_real_distribution<double> dist_zero_one(0.0,1.0);
  uniform_real_distribution<double> dist_azimut(0.0,2.0*Pi);
  uniform_real_distribution<double> dist_merging_section(0.0,interact_length);
  normal_distribution<double> dist_ion_beam_width(0.0,0.25*ion_beam_fwhm/sqrt(ln2)); // includes conversion from FWHM to Gaussian sigma

  int frag1_A_transmitted=0, frag2_A_transmitted=0, frag1_B_transmitted=0, frag2_B_transmitted=0;
  // *************************************  Monte-Carlo loop ************************
  for (int i=0; i<MonteCarlo_iterations;i++)
    {
      // pick a position in the x-y-plane of the beam 
      double b = dist_ion_beam_width(randgen);
      double phib = dist_azimut(randgen);      
      double bx=b*cos(phib);
      double by=b*sin(phib);

      // pick a random position on the beam axis in the interaction region
      double s = dist_merging_section(randgen);
      
      // pick a random molecular orientation
      double theta = random_theta(dist_zero_one(randgen),frag_anisotropy);
      double phi = dist_azimut(randgen);

      // just for checking against the analytical formula
      //double radiusA_frag1 = (drift_detA+s)*frag1_velocity_ratio*sin(theta);
      //fout << radiusA_frag1;    

      // fragment 1
      double s1r=s*frag1_velocity_ratio*sin(theta); 
      double s1x=s1r*cos(phi); 
      double s1y=s1r*sin(phi);
      double s1bx = s1x+bx;
      double s1by = s1y+by;
      double delta1 = frag1_velocity_ratio*cos(theta);
      double frag1_radius = sqrt(s1bx*s1bx+s1by*s1by);
      bool frag1_I_flag = ((frag1_radius < 0.5*dia_interact) && (fabs(s1bx) < 0.5*dia_trinos_x) && (fabs(s1by) < 0.5*dia_trinos_y));
      Vector5f v_frag1_int_region, v_frag1_detectorA, v_frag1_detectorB;
      v_frag1_int_region << s1bx, s1x/s, delta1, s1by, s1y/s;
      // check whether particle is on detector A
      v_frag1_detectorA = m_drift_detA*v_frag1_int_region;
      frag1_radius = sqrt(pow(v_frag1_detectorA(0),2)+pow(v_frag1_detectorA(3),2));
      bool frag1_A_flag = frag1_I_flag && (frag1_radius < 0.5*dia_detectorA);
      // check whether particle passes entrance aperture of demerging magnet
      v_frag1_detectorB = m_drift_demerger*v_frag1_int_region;
      frag1_radius = sqrt(pow(v_frag1_detectorB(0),2)+pow(v_frag1_detectorB(3),2));
      bool frag1_B_flag = frag1_I_flag && (frag1_radius < 0.5*dia_demerger);     
      
      // fragment 2
      double s2r=s*frag2_velocity_ratio*sin(Pi-theta); 
      double s2x=s2r*cos(Pi+phi); 
      double s2y=s2r*sin(Pi+phi);
      double s2bx = s2x+bx;
      double s2by = s2y+by;
      double delta2 = frag_momentum_ratio*cos(Pi-theta);
      double frag2_radius = sqrt(s2bx*s2bx+s2by*s2by);
      bool frag2_I_flag =  ((frag2_radius < 0.5*dia_interact) && (fabs(s2bx) < 0.5*dia_trinos_x) && (fabs(s2by) < 0.5*dia_trinos_y));
      Vector5f v_frag2_int_region, v_frag2_detectorA, v_frag2_detectorB;
      v_frag2_int_region << s2bx, s2x/s, delta2, s2by, s2y/s;
      // check whether particle is on detector A
      v_frag2_detectorA = m_drift_detA*v_frag2_int_region;
      frag2_radius = sqrt(pow(v_frag2_detectorA(0),2)+pow(v_frag2_detectorA(3),2));
      bool frag2_A_flag = frag2_I_flag && (frag2_radius < 0.5*dia_detectorA);
      // check whether particle passes entrance aperture of demerging magnet
      v_frag2_detectorB= m_drift_demerger*v_frag2_int_region;
      frag2_radius = sqrt(pow(v_frag2_detectorB(0),2)+pow(v_frag2_detectorB(3),2));
      bool frag2_B_flag = frag2_I_flag && (frag2_radius < 0.5*dia_demerger);

      // output to file
      if ((!frag1_A_flag) && (!frag2_A_flag) && (!frag1_B_flag) && (!frag2_B_flag)) continue; // no output to be written
      
      if (frag1_A_flag)
	{
	  frag1_A_transmitted++;
	  fout << v_frag1_detectorA(0) << ", " << v_frag1_detectorA(3)  << ", ";
	}
      else
	{
	  fout << " --, --,";
	}
      if (frag1_B_flag)
	{
	  frag1_B_transmitted++;
	  v_frag1_detectorB = m_drift_detB*m_demerger*m_drift_demerger*v_frag1_int_region;;
	  fout << v_frag1_detectorB(0) << ", " << v_frag1_detectorB(2)  << ", ";
	}
      else
	{
	  fout << " -- ,--, ";
	}
      if (frag2_A_flag)
	{
	  frag2_A_transmitted++;
	  fout << v_frag2_detectorA(0) << ", " << v_frag2_detectorA(3)  << ", ";
	}
      else
	{
	  fout << " --, --,";
	}
      if (frag2_B_flag)
	{
	  frag2_B_transmitted++;
	  v_frag2_detectorB = m_drift_detB*m_demerger*m_drift_demerger*v_frag2_int_region;;
	  fout << v_frag2_detectorB(0) << ", " << v_frag2_detectorB(2);
	}
      else
	{
	  fout << " -- ,--";
	}
      fout << endl;
    } // end for(i ...) end of Monte-Carlo loop
  //*******************************************
  fout << "#######################################################################################################################" << endl;
  fout << "###  fragm. 1 on detect. A : " << frag1_A_transmitted << endl;
  fout << "###  fragm. 1 on detect. B : " << frag1_B_transmitted << endl;
  fout << "###  fragm. 2 on detect. A : " << frag2_A_transmitted << endl;
  fout << "###  fragm. 2 on detect. B : " << frag2_B_transmitted << endl;
  fout.close();

  //if (!diagnostic_flag) { cout << "." << flush; }
  // ***************************** cleanup ******************************
  cout << endl << endl <<  " Data written to output file: " << outfilename << endl << endl;
  
  return 0;
}

