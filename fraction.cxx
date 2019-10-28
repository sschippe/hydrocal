/**
 * @file fraction.cxx
 *
 * @brief Surviving fractions of nl-Rydberg states populated in storage-ring recombination experiments
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: fraction.cxx 450 2017-08-29 08:38:26Z iamp $
 @endverbatim
 *
 */

#include <stdio.h>
#include <math.h>
#include <cstring>
#include <iostream>
#include "radrate.h"
#include "sigmarr.h"
#include "hydromath.h"
#include "fieldion.h"

using namespace std;

//#define printcasc 1

const double pi = 3.1415926536;
const double clight=2.99792458e10;
const double melectron=510999.6;
const double atomic_mass_unit=931494095.4; 

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Classical cut-off quantum number for a hydrogenic ion in an electric field
 * 
 * @param efield electric field strength in V/cm
 * @param z nuclear charge
 */
int nf(double efield, double z)
  { 
  return int(sqrt(sqrt(Fau*z*z*z/efield/9))+0.5);
  }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Mean cut-off quntum number for a hydrogenic ion in an electric field
 * 
 * Formula developed by Andreas Wolf, Habilitation Thesis, University of Heidelberg, 1992, page 75.
 *
 * @param tau ion-flight time (in s) from cooler to field-ionizing field 
 * @param z nuclear charge 
 * @param nf calssical cut-off quantum number (used as start value in iteration)
 */
int ngamma(double tau, double z, int nf)
   {
   double fac = pow(z,4)*tau*0.357*2.142e10;
   double n_new=nf, n_old=0;
   while (fabs(n_new-n_old)>1.E-4)
      {
      n_old = n_new;
      n_new = fac*exp(-1.36*log(n_old-1))*(0.481+0.68*log(n_old-1));
      n_new = exp(log(n_new)/3);
      }
   return int(n_new);
   }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at CSRm cooler
 *
 * Interactive generation of a hydrocal batch-file for the calculation of field-ionization 
 * survival probabilities behind the CSRm electron cooler
 */
void setup_batch_csrm(void)
   {
     // geometric quantities, all lengths are in cm
   const double rho =  760.0; // bending radius of dipole magnet;
   const double phi =   22.5; // bending angle of dipole magnet;
   const double z1_t = 400.0; // distance from beginning of cooler to toroid
   const double z2_t =  10.0; // distance from end of cooler to toroid
   const double dz_t =  90.0; // distance within toroid
   const double z_d =  800.0; // distance from toroid to dipole

   double amu,gamma,vel,ecool,bpar,efield,ef_t,zsum;
   int z,nmax,ncasc,softcut;
   char fn[200],fnroot[200], answer[2];
   FILE *fout;

   printf("\n");
   printf("\n Setup of a batch file for the calculation of nl-specific");
   printf("\n    detection probabilities at the CSRm electron cooler");
   printf("\n");
   printf("\n Give effective nuclear charge ...........: ");
   scanf("%d",&z);
   printf("\n Give atomic mass in amu .................: ");
   scanf("%lf",&amu);
   printf("\n Give cooling energy in eV ...............: ");
   scanf("%lf",&ecool);
   printf("\n Give magnetic guiding field in mT .......: ");
   scanf("%lf",&bpar);

   bpar *= 1E-7 ; // conversion of mT to Vs/(cm)^2
   gamma = 1.0+ecool/melectron;
   vel = clight*sqrt(1.0-1.0/(gamma*gamma)); // in cm/s
   ef_t = 0.5*bpar*vel; // E field in toroid n V/cm 

   efield = 1.036427E-12*gamma*vel*vel*amu/z/rho; // in Dipole

// calculate the length of the recombined ion's path through the dipole magnet
   double dz_d = rho*z/(z-1.0)*asin((z-1.0)/z*cos(phi*pi/180.0));
   
   printf("\n The relativistic gamma factor is: %9.8f\n",gamma);
   printf("\n The ion velocity is ............: %9.5g cm/s\n",vel);

   printf("\n The flight times are (in the ion's frame of reference)");
   
   printf("\n from the beginning of the cooler to the toroid : %6.2f ns",
          z1_t*1e9/vel/gamma);
   printf("\n from the end of the cooler to the toroid ......: %6.2f ns",
          z2_t*1e9/vel/gamma);
   printf("\n from the toroi to to the dipole ...............: %6.2f ns",
          z_d*1e9/vel/gamma);
   printf("\n through the cooler ............................: %6.2f ns",
          (z1_t-z2_t)*1e9/vel/gamma);
   printf("\n through the toroid ............................: %6.2f ns",
          dz_t*1e9/vel/gamma);
   printf("\n through the dipole ............................: %6.2f ns",
          dz_d*1e9/vel/gamma);
 
   double dt_t, dt_d;
   zsum = 0.5*(z1_t+z2_t);
   dt_t = zsum/vel/gamma;
   printf("\n from the cooler center to the toroid ..........: %6.2f ns",
          dt_t*1e9);
   zsum += z_d;
   dt_d = zsum/vel/gamma;
   printf("\n from the cooler center to the dipole ..........: %6.2f ns",
          dt_d*1e9);

   int nF_t  = nf(ef_t,z); 
   int nF_d  = nf(efield,z);
   printf("\n\n The motional electric fields and cut-off quantum numbers are ");
   printf("\n in the toroid ..................: %9.2f kV/cm",ef_t*1E-3);
   printf(", nF = %3d",nF_t);
   printf("\n in the charge analyzing dipole .: %9.2f kV/cm",efield*1E-3);
   printf(", nF = %3d",nF_d);


   printf("\n Give maximum main quantum number ........: ");
   scanf("%d",&nmax); 
   printf("\n Give number of cascade steps ............: ");
   scanf("%d",&ncasc);
   printf("\n Hard or soft cut-off ? (h/s) ............: ");
   scanf("%s",&answer);
   softcut = ((!strcmp(answer,"s")) || (!strcmp(answer,"S")));
   printf("\n Give filename for output (*.fnl) ........: ");
   scanf("%s",&fnroot);

   strcpy(fn,fnroot);
   strcat(fn,".bat");
   
   fout = fopen(fn,"w");   

   fprintf(fout,"4\n");

// field ionization data for the toroid

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f 100\n",z1_t,z2_t);
   fprintf(fout,"45\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_t);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_t);
   if (softcut) nF_t = 0;
   fprintf(fout,"%4d\n",nF_t);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"n\n");

   strcpy(fn,fnroot);
   strcat(fn,"_t");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");


// field ionization data for the charge analyzing dipole magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f %5.1f\n",z_d,z_d,rho);
   fprintf(fout,"%5.1f\n",phi);
   fprintf(fout,"n\n");
   fprintf(fout,"n\n");
   if (softcut) nF_d = 0;
   fprintf(fout,"%4d\n",nF_d);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

   fprintf(fout,"0\n");
   fprintf(fout,"0\n");
   fclose(fout);

   printf("\n\n run batch job by issuing the command: ");
   printf("nice hydrocal <%s.bat >%s.out &\n",fnroot,fnroot);
   printf("-------------------------------------------------");
   printf("-------------------------------------------------\n\n");

   }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at CSR cooler
 *
 * Interactive generation of a hydrocal batch-file for the calculation of field-ionization 
 * survival probabilities behind the CSR electron cooler
 */
void setup_batch_csr(void)
   {
     // geometric quantities, all lengths are in cm
   const double dz_c  =   100.0; // length of the cooler  
   const double z_c_to =   60.0; // distance from the center of the cooler to the toriod  
   const double dz_to  =    12.0; // length of orbit in the toriod
   const double toroid_angle = 0.4; // max deflection angle of the electron beam in the toroid in rad

   const double z_c_06  =  215.5; // distance from the center of the cooler to the 6-deg deflector  
   const double dz_06 =    20.3; // length of orbit in 06-deg deflector
   const double rho06 =   200.0; // bending radius of 6-deg deflector;
   const double gap06 =    12.0; // distance between bending plates of 6-deg deflector

   const double z_06_39 =  110.7; // distance between 6-deg and 39-deg deflector
   const double dz_39 =    67.7; // length of orbit in 39-deg deflector
   const double rho39 =   100.0; // bending radius of 39-deg deflector;
   const double gap39 =     6.0; // distance between bending plates of 39-deg deflector

   int zeff;
   double mion, eion, Bcooler;
   printf("\n");
   printf("\n Setup of a batch file for the calculation of nl-specific");
   printf("\n    detection probabilities at the CSR electron cooler");
   printf("\n");
   printf("\n Give effective nuclear charge ...........: ");
   scanf("%d",&zeff);
   printf("\n Give atomic mass in amu .................: ");
   scanf("%lf",&mion);
   printf("\n Give ion energy in keV ..................: ");
   scanf("%lf",&eion);
   printf("\n Give magnetic guiding field in mT .......: ");
   scanf("%lf",&Bcooler);
   Bcooler *= 1E-7 ; // conversion of mT to Vs/(cm)^2

   eion *= 1000.0;
   double ecool = eion*melectron/(mion*atomic_mass_unit);
   double gamma = 1.0+ecool/melectron;
   double vel = clight*sqrt(1.0-1.0/(gamma*gamma)); // in cm/s
   
   printf("\n The relativistic gamma factor is: %9.8f\n",gamma);
   printf("\n The ion velocity is ............: %9.5g cm/s\n",vel);
   printf("\n The cooling energy is ..........: %9.5g eV\n",ecool);

   printf("\n The flight times are (in the ion's frame of reference)");

   printf("\n through the cooler ............................: %6.2f mus", dz_c*1e6/vel/gamma);
   printf("\n through the toroid ............................: %6.2f mus", dz_to*1e6/vel/gamma);
   printf("\n through the  6-deg deflector ..................: %6.2f mus", dz_06*1e6/vel/gamma);
   printf("\n through the 39-deg deflector ..................: %6.2f mus", dz_39*1e6/vel/gamma);
 
   double dt_to = z_c_to/vel/gamma;
   printf("\n from the cooler center to the toroid ..........: %6.2f mus", dt_to*1e6);
   double dt_06 = z_c_06/vel/gamma;
   printf("\n from the cooler center to the  6-deg deflector.: %6.2f mus", dt_06*1e6);
   double dt_39 = (z_c_06+dz_06+z_06_39)/vel/gamma;
   printf("\n from the cooler center to the 39-deg deflector.: %6.2f mus", dt_39*1e6);

   double ef_to = Bcooler*sin(toroid_angle)*vel;
   double ef_06 = 2.0*eion/zeff*log((rho06+0.5*gap06)/(rho06-0.5*gap06));
   double ef_39 = 2.0*eion/zeff*log((rho39+0.5*gap39)/(rho39-0.5*gap39));
   int nF_to  = nf(ef_to,zeff);
   int nF_06 = nf(ef_06,zeff);
   int nF_39 = nf(ef_39,zeff);

   printf("\n\n The motional electric fields and cut-off quantum numbers are ");
   printf("\n in the cooler toriod    .................: %9.2f kV/cm",ef_to*1E-3);
   printf(", nF = %3d",nF_to);
   printf("\n in the  6-deg deflector .................: %9.2f kV/cm",ef_06*1E-3);
   printf(", nF = %3d",nF_06);
   printf("\n in the 39-deg deflector .................: %9.2f kV/cm",ef_39*1E-3);
   printf(", nF = %3d",nF_39);

   int nmax, ncasc;
   char fn[200],fnroot[200], answer[2];
   printf("\n Give maximum main quantum number ........: ");
   scanf("%d",&nmax); 
   printf("\n Give number of cascade steps ............: ");
   scanf("%d",&ncasc);
   printf("\n Hard or soft cut-off ? (h/s) ............: ");
   scanf("%s",&answer);
   int softcut = ((!strcmp(answer,"s")) || (!strcmp(answer,"S")));
   printf("\n Give filename for output (*.fnl) ........: ");
   scanf("%s",&fnroot);

   strcpy(fn,fnroot);
   strcat(fn,".bat");
   
   FILE *fout;
   fout = fopen(fn,"w");   

   fprintf(fout,"4\n");

// field ionization data for the toroid

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",zeff);
   fprintf(fout,"%6.2f\n",mion);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f 100\n", z_c_to+0.5*dz_c, z_c_to-0.5*dz_c);
   fprintf(fout,"45\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_to);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_to);
   if (softcut) nF_to = 0;
   fprintf(fout,"%4d\n",nF_to);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"n\n");

   strcpy(fn,fnroot);
   strcat(fn,"_t");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

// field ionization data for the 6-deg deflector

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",zeff);
   fprintf(fout,"%6.2f\n",mion);
   fprintf(fout,"%7.5f\n",ecool);
   fprintf(fout,"0\n");
   double z_to_06 = z_c_06-z_c_to; // distance from toroid to 6deg deflector
   fprintf(fout,"%5.1f %5.1f %5.1f\n",z_to_06, z_to_06,rho06);
   fprintf(fout,"-6.0\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_06);
   fprintf(fout,"%5.1f\n",gap06);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_06);
   if (softcut) nF_06 = 0;
   fprintf(fout,"%4d\n",nF_06);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   strcat(fn,"_06");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

// field ionization data for the 39-deg deflector

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",zeff);
   fprintf(fout,"%6.2f\n",mion);
   fprintf(fout,"%7.5f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f %5.1f\n",dz_06+z_06_39,dz_06+z_06_39,rho39);
   fprintf(fout,"-39.0\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_39);
   fprintf(fout,"%5.1f\n",gap39);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_39);
   if (softcut) nF_39 = 0;
   fprintf(fout,"%4d\n",nF_39);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   strcat(fn,"_39");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

   fprintf(fout,"0\n");
   fprintf(fout,"0\n");

   fclose(fout);

   printf("\n\n run batch job by issuing the command: ");
   printf("nice hydrocal <%s.bat >%s.out &\n",fnroot,fnroot);
   printf("-------------------------------------------------");
   printf("-------------------------------------------------\n\n");

   }


//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at TSR cooler
 *
 * Interactive generation of a hydrocal batch-file for the calculation of field-ionization 
 * survival probabilities behind the TSR electron cooler
 */
void setup_batch_tsrc(void)
   {
   const double rho =  115.0; // bending radius of dipole magnet;
   const double phi =   45.0; // bending angle of dipole magnet;
   const double z1_t = 170.0; // distance from beginning of cooler to toroid
   const double z2_t =  20.0; // distance from end of cooler to toroid
   const double dz_t =  50.0; // 56.5; // distance within toroid
   const double z_c1 =  82.5; // distance from toroid to KDX1
   const double dz_c1=  34.0; // distance within KDX1
   const double z_c2 =  50.5; // distance from KDX1 to KDX2
   const double dz_c2=  18.0; // distance within KDX2
   const double z_d =  244.0; // distance from KDX2 to dipole

   double amu,gamma,vel,ecool,bpar,efield,ef_t,ef_c1,ef_c2,zsum;
   int z,nmax,ncasc,softcut;
   char fn[200],fnroot[200], answer[2];
   FILE *fout;

   printf("\n");
   printf("\n Setup of a batch file for the calculation of nl-specific");
   printf("\n    detection probabilities at the TSR electron cooler");
   printf("\n");
   printf("\n Give effective nuclear charge ...........: ");
   scanf("%d",&z);
   printf("\n Give atomic mass in amu .................: ");
   scanf("%lf",&amu);
   printf("\n Give cooling energy in eV ...............: ");
   scanf("%lf",&ecool);
   printf("\n Give magnetic guiding field in mT .......: ");
   scanf("%lf",&bpar);

   bpar *= 1E-7 ; // conversion of mT to Vs/(cm)^2
   gamma = 1.0+ecool/melectron;
   vel = clight*sqrt(1.0-1.0/(gamma*gamma)); // in cm/s
   efield = bpar*vel; // in V/cm
  
   ef_t   = 0.5*efield; // 0.36*efield
   ef_c1  = 0.8*efield;
   ef_c2  = 1.6*efield;
   efield = 1.036427E-12*gamma*vel*vel*amu/z/rho;

// calculate the length of the recombined ion's path through the dipole magnet
   double dz_d = rho*z/(z-1.0)*asin((z-1.0)/z*cos(phi*pi/180.0));
   
   printf("\n The relativistic gamma factor is: %9.8f\n",gamma);
   printf("\n The ion velocity is ............: %9.5g cm/s\n",vel);

   printf("\n The flight times are (in the ion's frame of reference)");
   
   printf("\n from the beginning of the cooler to the toroid : %6.2f ns",
          z1_t*1e9/vel/gamma);
   printf("\n from the end of the cooler to the toroid ......: %6.2f ns",
          z2_t*1e9/vel/gamma);
   printf("\n from the toroid to the first correction magnet.: %6.2f ns",
          z_c1*1e9/vel/gamma);
   printf("\n from the first to the second correction magnet.: %6.2f ns",
          z_c2*1e9/vel/gamma);
   printf("\n from the second correction magnet to the dipole: %6.2f ns",
          z_d*1e9/vel/gamma);
   printf("\n through the cooler ............................: %6.2f ns",
          (z1_t-z2_t)*1e9/vel/gamma);
   printf("\n through the toroid ............................: %6.2f ns",
          dz_t*1e9/vel/gamma);
   printf("\n through the first correction magnet............: %6.2f ns",
          dz_c1*1e9/vel/gamma);
   printf("\n through the second correction magnet...........: %6.2f ns",
          dz_c2*1e9/vel/gamma);
   printf("\n through the dipole ............................: %6.2f ns",
          dz_d*1e9/vel/gamma);
 
   double dt_t, dt_c1, dt_c2, dt_d;
   zsum = 0.5*(z1_t+z2_t);
   dt_t = zsum/vel/gamma;
   printf("\n from the cooler center to the toroid ..........: %6.2f ns",
          dt_t*1e9);
   zsum += z_c1;
   dt_c1 = zsum/vel/gamma;
   printf("\n from the cooler center to the 1st correction m.: %6.2f ns",
          dt_c1*1e9);
   zsum += z_c2;
   dt_c2 = zsum/vel/gamma;
   printf("\n from the cooler center to the 2nd correction m.: %6.2f ns",
          dt_c2*1e9);
   zsum += z_d;
   dt_d = zsum/vel/gamma;
   printf("\n from the cooler center to the dipole ..........: %6.2f ns",
          dt_d*1e9);

   int nF_t  = nf(ef_t,z); 
   int nF_c1 = nf(ef_c1,z);
   int nF_c2 = nf(ef_c2,z);
   int nF_d  = nf(efield,z);
   printf("\n\n The motional electric fields and cut-off quantum numbers are ");
   printf("\n in the toroid ..................: %9.2f kV/cm",ef_t*1E-3);
   printf(", nF = %3d",nF_t);
   printf("\n in the first correction magnet .: %9.2f kV/cm",ef_c1*1E-3);
   printf(", nF = %3d",nF_c1);
   printf("\n in the second correction magnet : %9.2f kV/cm",ef_c2*1E-3);
   printf(", nF = %3d",nF_c2);
   printf("\n in the charge analyzing dipole .: %9.2f kV/cm",efield*1E-3);
   printf(", nF = %3d",nF_d);


   printf("\n Give maximum main quantum number ........: ");
   scanf("%d",&nmax); 
   printf("\n Give number of cascade steps ............: ");
   scanf("%d",&ncasc);
   printf("\n Hard or soft cut-off ? (h/s) ............: ");
   scanf("%s",&answer);
   softcut = ((!strcmp(answer,"s")) || (!strcmp(answer,"S")));
   printf("\n Give filename for output (*.fnl) ........: ");
   scanf("%s",&fnroot);

   strcpy(fn,fnroot);
   strcat(fn,".bat");
   
   fout = fopen(fn,"w");   

   fprintf(fout,"4\n");

// field ionization data for the toroid

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f 100\n",z1_t,z2_t);
   fprintf(fout,"45\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_t);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_t);
   if (softcut) nF_t = 0;
   fprintf(fout,"%4d\n",nF_t);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"n\n");

   strcpy(fn,fnroot);
   strcat(fn,"_t");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

// field ionization data for the first correction magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f 100\n",z_c1,z_c1);
   fprintf(fout,"45\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_c1);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_c1);
   if (softcut) nF_c1 = 0;
   fprintf(fout,"%4d\n",nF_c1);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   strcat(fn,"_c1");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

// field ionization data for the second correction magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f 100\n",z_c2,z_c2);
   fprintf(fout,"45\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_c2);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_c2);
   if (softcut) nF_c2 = 0;
   fprintf(fout,"%4d\n",nF_c2);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   strcat(fn,"_c2");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

// field ionization data for the charge analyzing dipole magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f %5.1f\n",z_d,z_d,rho);
   fprintf(fout,"%5.1f\n",phi);
   fprintf(fout,"n\n");
   fprintf(fout,"n\n");
   if (softcut) nF_d = 0;
   fprintf(fout,"%4d\n",nF_d);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

   fprintf(fout,"0\n");
   fprintf(fout,"0\n");
   fclose(fout);

   printf("\n\n run batch job by issuing the command: ");
   printf("nice hydrocal <%s.bat >%s.out &\n",fnroot,fnroot);
   printf("-------------------------------------------------");
   printf("-------------------------------------------------\n\n");

   }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at TSR target
 *
 * Interactive generation of a hydrocal batch-file for the calculation of field-ionization 
 * survival probabilities behind the TSR electron target
 */
void setup_batch_tsrt(void)
   {
   const double rho  =  115.0; // bending radius of dipole magnet;
   const double phi  =   45.0; // bending angle of dipole magnet;
   const double z1_t =  156.8; // distance from beginning of target to toroid (in this case 50 % of field)
   const double z2_t =   19.2; // end of target to beginning toroid (in this case 50 % of field)
   const double dz_t =   43.0; // distance within toroid (50% to 50%)
   const double z_c1 =   50.4; // distance from beginning toroid to beginning KDY2
                               //  (Yes, KDY2 is nearer to the toroide than KDY1) 
   const double dz_c1=   21.0; // distance within KDY2 (2cm each side linear fringe assumed so +1cm for each side)
   const double z_c2 =   35.5; // distance from beginning KDY2 to beginning KDY1 
   const double dz_c2=   12.0; // distance within KDY1 
   const double z_d  =  298.1; // distance from beginning KDY1 to beginning dipole

   double amu,gamma,vel,ecool,bpar,efield,ef_t,ef_c1,ef_c2,zsum;
   int z,nmax,ncasc,softcut;
   char fn[200],fnroot[200], answer[2];
   FILE *fout;

   printf("\n");
   printf("\n Setup of a batch file for the calculation of nl-specific");
   printf("\n    detection probabilities at the TSR electron target");
   printf("\n");
   printf("\n Give effective nuclear charge ...........: ");
   scanf("%d",&z);
   printf("\n Give atomic mass in amu .................: ");
   scanf("%lf",&amu);
   printf("\n Give cooling energy in eV ...............: ");
   scanf("%lf",&ecool);
   printf("\n Give magnetic guiding field in mT .......: ");
   scanf("%lf",&bpar);

   bpar *= 1E-7 ; // conversion of mT to Vs/(cm)^2
   gamma = 1.0+ecool/melectron;
   vel = clight*sqrt(1.0-1.0/(gamma*gamma)); // in cm/s
   efield = bpar*vel; // in V/cm
  
   ef_t  = 0.5*efield; // EWS
   ef_c1 = 1.9*efield; // EWS
   ef_c2 = 1.0*efield; // EWS
   efield = 1.036427E-12*gamma*vel*vel*amu/z/rho;

// calculate the length of the recombined ion's path through the dipole magnet
   double dz_d = rho*z/(z-1.0)*asin((z-1.0)/z*cos(phi*pi/180.0));
   
   printf("\n The relativistic gamma factor is: %9.8f\n",gamma);
   printf("\n The ion velocity is ............: %9.5g cm/s\n",vel);

   printf("\n The flight times are (in the ion's frame of reference)");
   
   printf("\n from the beginning of the target to the toroid : %6.2f ns",
          z1_t*1e9/vel/gamma);
   printf("\n from the end of the target to the toroid ......: %6.2f ns",
          z2_t*1e9/vel/gamma);
   printf("\n from the target to the first correction magnet.: %6.2f ns",
          z_c1*1e9/vel/gamma);
   printf("\n from the first to the second correction magnet.: %6.2f ns",
          z_c2*1e9/vel/gamma);
   printf("\n from the second correction magnet to the dipole: %6.2f ns",
          z_d*1e9/vel/gamma);
   printf("\n through the target ............................: %6.2f ns",
          (z1_t-z2_t)*1e9/vel/gamma);
   printf("\n through the toroid ............................: %6.2f ns",
          dz_t*1e9/vel/gamma);
   printf("\n through the first correction magnet............: %6.2f ns",
          dz_c1*1e9/vel/gamma);
   printf("\n through the second correction magnet...........: %6.2f ns",
          dz_c2*1e9/vel/gamma);
   printf("\n through the dipole ............................: %6.2f ns",
          dz_d*1e9/vel/gamma);
 
   double dt_t, dt_c1, dt_c2, dt_d;
   zsum = 0.5*(z1_t+z2_t);
   dt_t = zsum/vel/gamma;
   printf("\n from the target center to the toroid ..........: %6.2f ns",
          dt_t*1e9);
   zsum += z_c1;
   dt_c1 = zsum/vel/gamma;
   printf("\n from the target center to the 1st correction m.: %6.2f ns",
          dt_c1*1e9);
   zsum += z_c2;
   dt_c2 = zsum/vel/gamma;
   printf("\n from the target center to the 2nd correction m.: %6.2f ns",
          dt_c2*1e9);
   zsum += z_d;
   dt_d = zsum/vel/gamma;
   printf("\n from the target center to the dipole ..........: %6.2f ns",
          dt_d*1e9);


   int nF_t  = nf(ef_t,z); 
   int nF_c1 = nf(ef_c1,z);
   int nF_c2 = nf(ef_c2,z);
   int nF_d  = nf(efield,z);
   printf("\n\n The motional electric fields and cut-off quantum numbers are ");
   printf("\n in the toroid ..................: %9.2f kV/cm",ef_t*1E-3);
   printf(", nF = %3d",nF_t);
   printf("\n in the first correction magnet .: %9.2f kV/cm",ef_c1*1E-3);
   printf(", nF = %3d",nF_c1);
   printf("\n in the second correction magnet : %9.2f kV/cm",ef_c2*1E-3);
   printf(", nF = %3d",nF_c2);
   printf("\n in the charge analyzing dipole .: %9.2f kV/cm",efield*1E-3);
   printf(", nF = %3d",nF_d);


   printf("\n Give maximum main quantum number ........: ");
   scanf("%d",&nmax); 
   printf("\n Give number of cascade steps ............: ");
   scanf("%d",&ncasc);
   printf("\n Hard or soft cut-off ? (h/s) ............: ");
   scanf("%s",&answer);
   softcut = ((!strcmp(answer,"s")) || (!strcmp(answer,"S")));
   printf("\n Give filename for output (*.fnl) ........: ");
   scanf("%s",&fnroot);

   strcpy(fn,fnroot);
   strcat(fn,".bat");
   
   fout = fopen(fn,"w");   

   fprintf(fout,"4\n");

// field ionization data for the toroid

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f 100\n",z1_t,z2_t);
   fprintf(fout,"45\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_t);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_t);
   if (softcut) nF_t = 0;
   fprintf(fout,"%4d\n",nF_t);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"n\n");

   strcpy(fn,fnroot);
   strcat(fn,"_t");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

// field ionization data for the first correction magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f 100\n",z_c1,z_c1);
   fprintf(fout,"45\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_c1);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_c1);
   if (softcut) nF_c1 = 0;
   fprintf(fout,"%4d\n",nF_c1);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   strcat(fn,"_c1");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

// field ionization data for the second correction magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f 100\n",z_c2,z_c2);
   fprintf(fout,"45\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_c2);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_c2);
   if (softcut) nF_c2 = 0;
   fprintf(fout,"%4d\n",nF_c2);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   strcat(fn,"_c2");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

// field ionization data for the charge analyzing dipole magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f %5.1f\n",z_d,z_d,rho);
   fprintf(fout,"%5.1f\n",phi);
   fprintf(fout,"n\n");
   fprintf(fout,"n\n");
   if (softcut) nF_d = 0;
   fprintf(fout,"%4d\n",nF_d);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

   fprintf(fout,"0\n");
   fprintf(fout,"0\n");
   fclose(fout);

   printf("\n\n run batch job by issuing the command: ");
   printf("nice hydrocal <%s.bat >%s.out &\n",fnroot,fnroot);
   printf("-------------------------------------------------");
   printf("-------------------------------------------------\n\n");

   }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at ESR cooler
 *
 * Interactive generation of a hydrocal batch-file for the calculation of field-ionization 
 * survival probabilities behind the ESR electron cooler
 */
void setup_batch_esr(void)
{
   const double rho =    625.0; // bending radius of dipole magnet;
   const double phi =     60.0; // bending angle of dipole magnet;
   const double z1_t =   295.0; // distance from beginning of cooler to toroid
   const double z2_t =    45.0; // distance from end of cooler to toroid
   const double dz_t =    60.0; // distance within toroid
   const double z_d1 =   768.0; // distance from toroid to 1st dipole
   const double z_d2 = 1443.75; // distance from 1st dipole to 2nd dipole
   const double z_d3 = 1443.75; // distance from 2st dipole to 3rd dipole
   const double z_d4 = 2530.50; // distance from 3rd dipole to 4th dipole
   const double dz_d =  654.50; // distance within each dipole


   double amu,gamma,vel,ecool,bpar,efield,ef_t,zsum,z1,z2;
   int z,nmax,ncasc,ndip,softcut,use_toroid;
   char fn[200],fnroot[200],answer[2];
   FILE *fout;

   printf("\n");
   printf("\n Setup of a batch file for the calculation of nl-specific");
   printf("\n    detection probabilities at the ESR electron cooler");
   printf("\n");
   printf("\n Give effective nuclear charge ................: ");
   scanf("%d",&z);
   printf("\n Give atomic mass in amu ......................: ");
   scanf("%lf",&amu);
   printf("\n Give cooling energy in eV ....................: ");
   scanf("%lf",&ecool);
   printf("\n Give magnetic guiding field in mT ............: ");
   scanf("%lf",&bpar);
   printf("\n Give number of bending dipoles before detector: ");
   scanf("%d",&ndip);

   bpar *= 1E-7 ; // conversion of mT to Vs/(cm)^2
   gamma = 1.0+ecool/melectron;
   vel = clight*sqrt(1.0-1.0/pow(gamma,2)); // in cm/s
   efield = bpar*vel; // in V/cm
  
   ef_t   = 0.5*efield; // 0.36*efield
   efield = 1.036427E-12*gamma*vel*vel*amu/z/rho;

// calculate the length of the recombined ion's path through the dipole magnet
// double  dz_d = rho*z/(z-1.0)*asin((z-1.0)/z*cos(phi*pi/180.0));
   
   printf("\n The relativistic gamma factor is: %9.8f\n",gamma);
   printf("\n The ion velocity is ............: %9.5g cm/s\n",vel);

   printf("\n The flight times are (in the ion's frame of reference)");
   
   printf("\n from the beginning of the cooler to the toroid : %6.2f ns",
          z1_t*1e9/vel/gamma);
   printf("\n from the end of the cooler to the toroid ......: %6.2f ns",
          z2_t*1e9/vel/gamma);
   printf("\n from the toroid to the 1st dipole .............: %6.2f ns",
          z_d1*1e9/vel/gamma);
   if (ndip>1)
     {
     printf("\n from the 1st to the 2nd dipole ................: %6.2f ns",
            z_d2*1e9/vel/gamma);
     }
   if (ndip>2)
     {
     printf("\n from the 2nd to the 3rd dipole ................: %6.2f ns",
            z_d3*1e9/vel/gamma);
     }
   if (ndip>3)
     {
     printf("\n from the 3rd to the 4th dipole ................: %6.2f ns",
            z_d4*1e9/vel/gamma);
     }
   printf("\n through the cooler ............................: %6.2f ns",
          (z1_t-z2_t)*1e9/vel/gamma);
   printf("\n through the toroid ............................: %6.2f ns",
          dz_t*1e9/vel/gamma);
   printf("\n through the dipole ............................: %6.2f ns",
          dz_d*1e9/vel/gamma);



   double dt_t, dt_d1, dt_d2, dt_d3, dt_d4;
   zsum = 0.5*(z1_t+z2_t);
   dt_t = zsum/vel/gamma;
   printf("\n from the cooler center to the toroid ..........: %6.2f ns",
          dt_t*1e9);
   zsum += z_d1;
   dt_d1 = zsum/vel/gamma;
   printf("\n from the cooler center to the 1st dipole ......: %6.2f ns",
          dt_d1*1e9);
   if (ndip>1)
     {
     zsum += z_d2;
     dt_d2 = zsum/vel/gamma;
     printf("\n from the cooler center to the 2nd dipole ......: %6.2f ns",
            dt_d2*1e9);
     }
   if (ndip>2)
     {
     zsum += z_d3;
     dt_d3 = zsum/vel/gamma;
     printf("\n from the cooler center to the 3rd dipole ......: %6.2f ns",
            dt_d3*1e9);
     }
   if (ndip>3)
     {
     zsum += z_d4;
     dt_d4 = zsum/vel/gamma;
     printf("\n from the cooler center to the 4th dipole ......: %6.2f ns",
            dt_d4*1e9);
     }

   int nF_t = nf(ef_t,z);  
   int nF_d = nf(efield,z);
   printf("\n\n The motional electric fields and cut-off quantum numbers are ");
   printf("\n in the toroid ..................: %9.2f kV/cm",ef_t*1E-3);
   printf(", nF = %3d",nF_t);
   printf("\n in the charge analyzing dipole .: %9.2f kV/cm",efield*1E-3);
   printf(", nF = %3d",nF_d);


   printf("\n\n Consider toroid ? (y/n) .................: ");
   scanf("%s",&answer);
   use_toroid = ((!strcmp(answer,"y")) || (!strcmp(answer,"Y")));

   printf("\n Give maximum main quantum number ........: ");
   scanf("%d",&nmax); 
   printf("\n Give number of cascade steps ............: ");
   scanf("%d",&ncasc);
   printf("\n Hard or soft cut-off ? (h/s) ............: ");
   scanf("%s",&answer);
   softcut = ((!strcmp(answer,"s")) || (!strcmp(answer,"S")));
   printf("\n Give filename for output (*.fnl) ........: ");
   scanf("%s",&fnroot);

   strcpy(fn,fnroot);
   strcat(fn,".bat");
   
   fout = fopen(fn,"w");   

   fprintf(fout,"4\n");

// field ionization data for the toroid

   if (use_toroid)
     {
     fprintf(fout,"3\n");
     fprintf(fout,"%2d\n",z);
     fprintf(fout,"%6.2f\n",amu);
     fprintf(fout,"%9.2f\n",ecool);
     fprintf(fout,"0\n");
     fprintf(fout,"%5.1f %5.1f 100\n",z1_t,z2_t);
     fprintf(fout,"45\n");
     fprintf(fout,"y\n");
     fprintf(fout,"%5.1f\n",dz_t);
     fprintf(fout,"y\n");
     fprintf(fout,"%9.2f\n",ef_t);
     if (softcut) nF_t = 0;
     fprintf(fout,"%4d\n",nF_t);
     fprintf(fout,"%4d\n",nmax);
     fprintf(fout,"y\n");
     fprintf(fout,"%4d\n",ncasc);
     fprintf(fout,"n\n");

     strcpy(fn,fnroot);
     strcat(fn,"_t");
     fprintf(fout,"%s\n",fn);
     fprintf(fout,"n\n");
     z1 = 0.0; z2 = 0.0;
     }
   else
     {
     z1 = z1_t; z2= z2_t;
     }

// field ionization data for the 1st dipole magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f %5.1f\n",z_d1+z1,z_d1+z2,rho);
   fprintf(fout,"%5.1f\n",phi);
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_d);
   fprintf(fout,"n\n");
   if (softcut) nF_d = 0;
   fprintf(fout,"%4d\n",nF_d);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   if (use_toroid)
     {
     fprintf(fout,"y\n");
     fprintf(fout,"%s\n",fn);
     }
   else
     {
     fprintf(fout,"n\n");
     }

   strcpy(fn,fnroot);
   if (ndip>1) strcat(fn,"_d1");

   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

   if (ndip<2) goto finish;

// field ionization data for the 2nd dipole magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f %5.1f\n",z_d2,z_d2,rho);
   fprintf(fout,"%5.1f\n",phi);
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_d);
   fprintf(fout,"n\n");
   if (softcut) nF_d = 0;
   fprintf(fout,"%4d\n",nF_d);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   if (ndip>2) strcat(fn,"_d2");

   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

   if (ndip<3) goto finish;

// field ionization data for the 3rd dipole magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f %5.1f\n",z_d3,z_d3,rho);
   fprintf(fout,"%5.1f\n",phi);
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_d);
   fprintf(fout,"n\n");
   if (softcut) nF_d = 0;
   fprintf(fout,"%4d\n",nF_d);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   if (ndip>3) strcat(fn,"_d3");

   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

   if (ndip<4) goto finish;

// field ionization data for the 4th dipole magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f %5.1f\n",z_d4,z_d4,rho);
   fprintf(fout,"%5.1f\n",phi);
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_d);
   fprintf(fout,"n\n");
   if (softcut) nF_d = 0;
   fprintf(fout,"%4d\n",nF_d);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);

   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");


finish:

   fprintf(fout,"0\n");
   fprintf(fout,"0\n");

   fclose(fout);

   printf("\n\n run batch job by issuing the command: ");
   printf("nice hydrocal <%s.bat >%s.out &\n",fnroot,fnroot);
   printf("-------------------------------------------------");
   printf("-------------------------------------------------\n\n");
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief field ionization at CRYRING cooler
 *
 * Interactive generation of a hydrocal batch-file for the calculation of field-ionization 
 * survival probabilities behind the CRYRING electron cooler
 */
void setup_batch_cryring(void)
   {
   const double rho =  120.0; // bending radius of dipole magnet;
   const double phi =   30.0; // bending angle of dipole magnet;
   const double z1_t = 131.0; // distance from beginning of cooler to toroid
   const double z2_t =  27.0; // distance from end of cooler to toroid
   const double dz_t =  39.0; // distance within toroid
   const double z_c1 =  59.0; // distance from toroid to correction magent
   const double dz_c1=  21.0; // distance within correction magnet
   const double z_d =   47.0; // distance from correction magnet to dipole

   double amu,gamma,vel,ecool,bpar,efield,ef_t,ef_c1,ef_c2,zsum;
   int z,nmax,ncasc,softcut;
   char fn[200],fnroot[200],answer[2];
   FILE *fout;

   printf("\n");
   printf("\n Setup of a batch file for the calculation of nl-specific");
   printf("\n    detection probabilities at the cryring electron cooler");
   printf("\n");
   printf("\n Give effective nuclear charge ...........: ");
   scanf("%d",&z);
   printf("\n Give atomic mass in amu .................: ");
   scanf("%lf",&amu);
   printf("\n Give cooling energy in eV ...............: ");
   scanf("%lf",&ecool);
   printf("\n Give magnetic guiding field in mT .......: ");
   scanf("%lf",&bpar);

   bpar *= 1E-7 ; // conversion of mT to Vs/(cm)^2
   gamma = 1.0+ecool/melectron;
   vel = clight*sqrt(1.0-1.0/(gamma*gamma)); // in cm/s
   efield = bpar*vel; // in V/cm
  
   ef_t   = 0.5*efield;
   ef_c1  = 1.415*efield;
   efield = 1.036427E-12*gamma*vel*vel*amu/z/rho;

// calculate the length of the recombined ion's path through the dipole magnet
   double  dz_d = rho*z/(z-1.0)*asin((z-1.0)/z*cos(phi*pi/180.0));
   
   printf("\n The relativistic gamma factor is: %9.8f\n",gamma);
   printf("\n The ion velocity is ............: %9.5g cm/s\n",vel);

   printf("\n The flight times are (in the ion's frame of reference)");
   
   printf("\n from the beginning of the cooler to the toroid : %6.2f ns",
          z1_t*1e9/vel/gamma);
   printf("\n from the end of the cooler to the toroid ......: %6.2f ns",
          z2_t*1e9/vel/gamma);
   printf("\n from the toroid to the correction magnet ......: %6.2f ns",
          z_c1*1e9/vel/gamma);
   printf("\n from the correction magnet to the dipole ......: %6.2f ns",
          z_d*1e9/vel/gamma);
   printf("\n through the cooler ............................: %6.2f ns",
          (z1_t-z2_t)*1e9/vel/gamma);
   printf("\n through the toroid ............................: %6.2f ns",
          dz_t*1e9/vel/gamma);
   printf("\n through the correction magnet .................: %6.2f ns",
          dz_c1*1e9/vel/gamma);
   printf("\n through the dipole ............................: %6.2f ns",
          dz_d*1e9/vel/gamma);

   double dt_t, dt_c1, dt_d;
   zsum = 0.5*(z1_t+z2_t);
   dt_t = zsum/vel/gamma;
   printf("\n from the cooler center to the toroid ..........: %6.2f ns",
          dt_t*1e9);
   zsum += z_c1;
   dt_c1 = zsum/vel/gamma;
   printf("\n from the cooler center to the correction m.....: %6.2f ns",
          dt_c1*1e9);
   zsum += z_d;
   dt_d = zsum/vel/gamma;
   printf("\n from the cooler center to the dipole ..........: %6.2f ns",
          dt_d*1e9);

   int nF_t  = nf(ef_t,z); 
   int nF_c1 = nf(ef_c1,z);
   int nF_d  = nf(efield,z);
   printf("\n\n The motional electric fields are ");
   printf("\n in the toroid ..................: %9.2f kV/cm",ef_t*1E-3);
   printf(", nF = %3d",nF_t);
   printf("\n in the correction magnet .......: %9.2f kV/cm",ef_c1*1E-3);
   printf(", nF = %3d",nF_c1);
   printf("\n in the charge analyzing dipole .: %9.2f kV/cm",efield*1E-3);
   printf(", nF = %3d",nF_d);

   printf("\n\n Give maximum main quantum number ........: ");
   scanf("%d",&nmax); 
   printf("\n Give number of cascade steps ............: ");
   scanf("%d",&ncasc);
   printf("\n Hard or soft cut-off ? (h/s) ............: ");
   scanf("%s",&answer);
   softcut = ((!strcmp(answer,"s")) || (!strcmp(answer,"S")));
   printf("\n Give filename for output (*.fnl) ........: ");
   scanf("%s",&fnroot);

   strcpy(fn,fnroot);
   strcat(fn,".bat");
   
   fout = fopen(fn,"w");   

   fprintf(fout,"4\n");

// field ionization data for the toroid

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f 100\n",z1_t,z2_t);
   fprintf(fout,"45\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_t);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_t);
   if (softcut) nF_t = 0;
   fprintf(fout,"%4d\n",nF_t);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"n\n");

   strcpy(fn,fnroot);
   strcat(fn,"_t");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

// field ionization data for the correction magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f 100\n",z_c1,z_c1);
   fprintf(fout,"45\n");
   fprintf(fout,"y\n");
   fprintf(fout,"%5.1f\n",dz_c1);
   fprintf(fout,"y\n");
   fprintf(fout,"%9.2f\n",ef_c1);
   if (softcut) nF_c1 = 0;
   fprintf(fout,"%4d\n",nF_c1);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   strcat(fn,"_c");
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

// field ionization data for the charge analyzing dipole magnet

   fprintf(fout,"3\n");
   fprintf(fout,"%2d\n",z);
   fprintf(fout,"%6.2f\n",amu);
   fprintf(fout,"%9.2f\n",ecool);
   fprintf(fout,"0\n");
   fprintf(fout,"%5.1f %5.1f %5.1f\n",z_d,z_d,rho);
   fprintf(fout,"%5.1f\n",phi);
   fprintf(fout,"n\n");
   fprintf(fout,"n\n");
   if (softcut) nF_d = 0;
   fprintf(fout,"%4d\n",nF_d);
   fprintf(fout,"%4d\n",nmax);
   fprintf(fout,"y\n");
   fprintf(fout,"%4d\n",ncasc);
   fprintf(fout,"y\n");
   fprintf(fout,"%s\n",fn);

   strcpy(fn,fnroot);
   fprintf(fout,"%s\n",fn);
   fprintf(fout,"n\n");

   fprintf(fout,"0\n");
   fprintf(fout,"0\n");
   fclose(fout);

   printf("\n\n run batch job by issuing the command: ");
   printf("nice hydrocal <%s.bat >%s.out &\n",fnroot,fnroot);
   printf("-------------------------------------------------");
   printf("-------------------------------------------------\n\n");

   }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Entry to interactive batch-file setup for field-ionization calculations
 */
void setup_batch(void)
{
  int choice;

  start:
  printf("\n Which storage ring (target device)?\n");
  printf("\n 1) TSR cooler");
  printf("\n 2) TSR target");
  printf("\n 3) ESR");
  printf("\n 4) CRYRING");
  printf("\n 5) CSR");
  printf("\n 6) CSRm\n");
  printf("\n Make your choice : ");
  scanf("%d",&choice);

  switch(choice) {
    case 1: setup_batch_tsrc();break;
    case 2: setup_batch_tsrt();break;
    case 3: setup_batch_esr();break;
    case 4: setup_batch_cryring();break;
    case 5: setup_batch_csr();break;
    case 6: setup_batch_csrm();break;
    default: goto start;
  }
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Reads previously calculated survival probabilities from file
 */
int readfraction(double* fraction, char *header, int nmax)
  {
  FILE* fin;
  char fn[200], dummy[200];

  do
    {
    printf("\n Give filename for fraction table (*.fnl) ..: ");
    scanf("%s",&fn);
    strcat(fn,".fnl");
    fin = fopen(fn,"r");
    } while (!fin);
  fgets(header,200,fin);
  fgets(dummy,200,fin);
  fgets(dummy,200,fin);
  int n=0,l;
  double val1, val2 ,val3;
  while (n<=nmax)
     {
     if (fscanf(fin,"%d %d %lf %lf %lf",&n,&l,&val1,&val2,&val3)==EOF) {break;}
     if (n>nmax) break;
     fraction[(n-1)*n/2+l] *= val1;
     }
  fclose(fin); 
  return n >= nmax ? nmax : n-1;
  }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Interactive calculation of field-ionization cut-off quantum numbers
 */ 
void ngamma(void)
   {
   const double clight=2.99792458e10;
   const double melectron=510999.6;
   double z, mass, ecool, cool_dip, rho, gamma, vel, efield;
   int nF, selection;
   char answer[2];

   printf("\n Calculation of maximum main quantum numbers which can");
   printf("\n contribute to RR in a storage ring measurement.\n");
   printf("\n The cut off due to field ionization in the dipole magnet");
   printf("\n is calulated as  nf = [z^3/(9F)]^1/4 where F is the");
   printf("\n motional electric field in V/cm.\n");
   printf("\n Also calculated is ngamma, i.e. the maximum main quantum");
   printf("\n number of Rydberg levels which on the way from the cooler");
   printf("\n to the dipole magnet radiatively decay to below nf and");
   printf("\n therefore can contribute to RR. It is calculated by");
   printf("\n iteratively solving the equation");
   printf("\n ngamma = {0.357*z^4*g0*t*[0.481 + ln(ngamma-1)]");
   printf("*(ngamma-1)^-1.38}^1/3");
   printf("\n where g0 = 2.142e10/s, and z and t");
   printf("\n denote the effective ionic charge and the flight time");
   printf("\n between cooler and magnet, respectively");
   printf("\n (see Habiltation thesis of A. Wolf).\n\n\n"); 
   printf("\n Give effective ionic charge z .................: ");
   scanf("%lf",&z);
   printf("\n Give atomic mass (amu) ........................: ");
   scanf("%lf",&mass);
   printf("\n Give cooling energy (eV) ......................: ");
   scanf("%lf",&ecool);
   printf("\n Give the distance D between the cooler center and");
   printf("\n entrance of the dipole magnet and the bending");
   printf("\n radius R of the magnet\n"); 
   printf("                         manual input : 0\n");
   printf("     TSR values (D= 472 cm, R=115 cm) : 1\n");
   printf("     ESR values (D= 938 cm, R=625 cm) : 2\n");
   printf(" CRYRING values (D= 185 cm, R=120 cm) : 3\n");
   printf("    CSRm values (D=1100 cm, R=760 cm) : 4\n");
   printf("    CSRe values (D=1086 cm, R=600 cm) : 5\n\n");
   printf(" Make a selection ....................: ");
   scanf("%d",&selection);
   switch(selection)
      {
      case 1: cool_dip =  472; rho=115; break;
      case 2: cool_dip =  938; rho=625; break;
      case 3: cool_dip =  185; rho=120; break;
      case 4: cool_dip = 1100; rho=760; break;
      case 5: cool_dip = 1086; rho=600; break;
      default:
         printf("\n Give distance D and radius R in cm: ");
         scanf("%lf %lf",&cool_dip,&rho);
      }
   gamma = 1.0 +ecool/melectron;
   vel = clight*sqrt(1.0-1.0/(gamma*gamma));
   printf("\n The ion velocity is                          v = %10.4g cm/s",
          vel);
   double tau = cool_dip/vel/gamma;
   printf("\n The flight time between cooler and magnet is t = %10.4g s\n",
          tau);
   efield = 1.036427E-12*vel*vel*mass/z/rho;
   for (;;)
     {
     nF = nf(efield,z);
     printf("\n The electric field in the dipole magnet is   F = %10.4g V/cm",
            efield);
     printf("\n leading to a cutoff main quantum number     nF = %10d",nF);
     printf("\n\n Is the electric field ok? (y/n) .....: ");
     scanf("%s",&answer); 
     if (!strcmp(answer,"y") || !strcmp(answer,"Y") ) {break;}
     printf("\n Give electric field in V/cm .........: ");
     scanf("%lf",&efield);
     }

   int ng = ngamma(tau,z,nF);
   if (ng>nF)
     {
     printf("\n higher Ryberg levels up to              ngamma = %10d",ng);
     printf("\n can decay radiatively to below              nF = %10d",nF); 
     }
   else
     {
     printf("\n ATTENTION: ngamma = %5d <= nF = %5d",ng,nF);
     printf("\n Use nmax = nF-1 = %5d as maximum quantum number",nF-1);
     printf(" contributing to RR.");
     }
   }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief recursive calculation of cascade contributions to field-ionization survival probabilities
 */
double cascade(int counter, int ncascstep, int n1, int l1, int n_zero,
               int *n_list, int *l_list, double *pd_list, double *lt_list,
               const RADRATE& hydro, double *psurv, double *pcascdecay)
  {
  counter++;

  n_list[counter] = n1;
  l_list[counter] = l1;
  pd_list[counter]  = hydro.pdecay(n1,l1);
  lt_list[counter]  = hydro.life(n1,l1);

  static int i, k, n2max, docascades;
  static double prod, fcascade;
  
  double pdecay = 0.0;
  for (i=0; i<=counter; i++)
    {
    prod = 1.0;
    for (int k=0; k<=counter; k++) 
      { 
        if (k!=i) prod *= 1.0-lt_list[k]/lt_list[i];
      } 
    pdecay += pd_list[i]/prod;
    }
 
  *pcascdecay = pdecay; // fraction of flux going into cascades

  double fnl = 0.0;
  if (counter<ncascstep) { n2max = n1-1;} // full calc. for remaining casc.
             else { n2max = (n1<=n_zero) ? n1-1 : n_zero;}
  for(int n2=n2max; n2>0; n2--)
    {
    for(int l2=l1-1; l2<=l1+1; l2+=2)
      {
      if ( ( l2<0) || (l2>=n2) ) continue;
      double pdcasc = 0.0;
      fcascade = 0.0;
      if (counter<ncascstep)
         {
         fcascade = cascade(counter, ncascstep, n2, l2, n_zero, n_list, l_list,
                           pd_list, lt_list, hydro, psurv, &pdcasc);
         }
      fnl += hydro.branch(n1,l1,n2,l2)*((pdecay-pdcasc)*psurv[n2*(n2-1)/2+l2]+fcascade);
      } // end for (l2...)
    } // end (for (n2...)

#ifdef printcasc
  for(i=0;i<counter;i++)
    {
    printf("%3d %3d -> ",n_list[i],l_list[i]);
    }
  printf(" %3d %3d : %8.5f\n",n1,l1,fnl);
#endif

  return fnl;
  }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Calculation of nl-dependent hydrogenic field-ionization survival probabilities
 * 
 *
 * @param iselect = 0: elaborate treatment iusing the formulas given in the appendix
 *                     of Schippers et al. ApJ 555, 1027 (2001)
 *
 * @param iselect = 1: simple method of Zong et al. JPB 31, 3729 (1998) 
 *
 * Several output files will be generated.
 * The output files are *.fn and *.fnl contain l-averaged and l-resolved survival 
 *    probabilities, respectively.
 * The output file *.fm1 contains l-resolved survival probabilities in matrix form for 
 *    import into graphics software.
 * The output file *.fm2 contains l-resolved survival probabilities weighted with the corresponing
 *    RR cross section in matrix form for import into graphics software.
 * The output file *.fm3 contains l-resolved hydrogenic decay probabilities in matrix form for 
 *    import into graphics software.
 *
 * See appendix of Schippers et al. ApJ 555, 1027 (2001) for further details.
 */

void fractionTable(int iselect)
  {
  const double clight=2.99792458e10;
  const double melectron=510999.6;
  const double pi = 3.1415926536;
  double z,x1,x2,xd,ecool,rho,phi,mass;
  int nselective, ncut, nmax, softcut, ncascstep, selection; 
  char answer[2], filenameroot[200], filename[200], pfn[200];
  FILE *fin=NULL, *fout=NULL, *fout2=NULL, *fmatrix1=NULL, *fmatrix2=NULL, *fmatrix3=NULL;

if (iselect)
  {
  printf("Calculation of detection probabilities of hydrogenic n,l Rydberg\n");
  printf("states in storage ring recombination experiments taking into\n");
  printf("account radiative decay on the way from the electron cooler to\n");
  printf("the dipole magnet.\n");
  printf("The simple formula of Zong et al., JPB 31, 3729 (1998) is used.\n");
  }
else
 {
  printf("Calculation of detection probabilities of hydrogenic n,l Rydberg\n");
  printf("states in storage ring recombination experiments taking into\n");
  printf("account radiative decay on the way from the electron cooler to\n");
  printf("the dipole magnet or electrostatic deflector where the (motional)\n");
  printf("electric field determines the survival probability.\n");
  printf("Exact hydrogenic transition rates and semi-empirical hydrogenic\n");
  printf("field ionization  rates are used.\n"); 
  printf("Cascading can be taken into account.\n\n");
  }
  printf("The fractions are written to a file (*.fnl) that can be used \n");
  printf("for the calculation of recombination cross sections or rates.\n"); 
  printf("l-averaged fractions are written to a separate file (*.fn).\n");
  printf("The n,l-selective output file *.fnl consists of 5 columns.\n");
  printf(" col 1: main quantum number n\n");
  printf(" col 2: angular momentum quantum number l\n");
  printf(" col 3: fraction of n,l-population contributing to recombination\n");
  printf(" col 4: QM cross section for RR into n,l shell at 1e-12 eV\n");
  printf(" col 5: col3*col4/(maximum cross section per n)\n\n");
  printf("Three further outputfiles *.fm1, *.fm2 and  *fm3 contain the\n"); 
  printf(" cols 3 and 5 and decay probabilities in matrix form,\n");
  printf("respectively. These can be used to generate 3D-plots\n");
  printf("e.g. by import into a matrix within Origin.\n\n");

  printf("\n Give effective nuclear charge ...................: ");
  scanf("%lf",&z);
  printf("\n Give atomic mass (amu) ..........................: ");
  scanf("%lf",&mass);
  printf("\n Give cooling energy (eV) ........................: ");
  scanf("%lf",&ecool);
  double gamma = 1.0+ecool/melectron;
  double vel = clight*sqrt(1.0-1.0/(gamma*gamma));
  double eion = ecool*mass*atomic_mass_unit/melectron;
  printf("\n The ion energy   is %12.5g keV.",eion*0.001);
  printf("\n The ion velocity is %12.5g cm/s.\n\n",vel);
  printf(" Input of geometry values");
  printf("\n  D1: distance between beginning of cooler");
  printf(" and dipole magnet / deflector entrance");
  printf("\n  D2: distance between end of cooler");
  printf(" and dipole magnet / deflector entrance");
  printf("\n   R: bending radius of dipole magnet / deflector");
  printf("\n PHI: deflection angle of dipole magnet / deflector\n\n");
  printf("   (D1, D2, and R in cm, PHI in deg) manual input : 0\n");
  printf("     TSR values : ................................: 1\n");
  printf("     ESR values  .................................: 2\n");
  printf(" CRYRING values  .................................: 3\n");
  printf("     CSR values (no dipole magnet) ...............: 4\n");
  printf("    CSRm values .................. ...............: 5\n");
  printf("    CSRe values .................. ...............: 6\n");
  printf(" Give selection ..................................: ");
  scanf("%d",&selection);
  switch(selection)
     {
     case 1: // TSR
        x1  = 547;
        x2  = 397;
        rho = 115;
        phi = 45;
        break;
     case 2: // ESR
        x1  = 1063;
        x2  =  813;
        rho =  625;
        phi = 60;
        break;        
     case 3: // CRYRING
        x1  =  237;
        x2  =  133;
        rho =  120;
        phi = 30;
        break;        
     case 4: //CSR Heidelberg
        x1  =   200;
        x2  =   100;
        rho =   100;
        phi = -39.0; // negative value indicates electrostatic deflection
        break;        
     case 5:  //CSRm Lanzhou [Xi et al., NIMA 488 (2012) 11]
        x1  = 1300;
        x2  =  900;
        rho =  760;
        phi = 22.5;
        break;        
     case 6:  //CSRe Lanzhou [Xi et al., NIMA 488 (2012) 11]
        x1  = 1286;
        x2  =  886;
        rho =  600;
        phi = 22.5;
        break;        
     default:
        printf("\n Give D1, D2, and R in cm ......................: ");
        scanf("%lf %lf %lf",&x1,&x2,&rho);
        printf("\n Give Phi in deg (PHI<0 for el.-stat. defl.) ...: ");
        scanf("%lf",&phi);
     }
  if (iselect==1) 
  {   // use simple formula of Zong et al.
      double xmean = 0.5*(x1+x2);
      x1=xmean;
      x2=xmean;
      printf("\n Using D = %5.1f cm for the simple formula of Zong et al.\n",xmean); 
  }

// calculate the length of the recombined ion's path through the dipole magnet
  if (phi>0)
    { //  magnetic deflection
      if (z==1)
	{
	  xd = 80;
	}
      else
	{
	  xd = rho*z/(z-1.0)*asin((z-1.0)/z*cos(phi*pi/180.0));
	}
  printf("\n The distance within the magnet is %5.3f cm.\n",xd);    }
  else
    { // electrostatic deleflection
      xd = fabs(rho*phi*pi/180.0);
      printf("\n The distance within the deflector is %5.3f cm.\n",xd);
    }

  printf("\n Do you want to change this value ? (y/n) ........: ");
  scanf("%s",answer);
  if ((!strcmp(answer,"y"))||(!strcmp(answer,"Y")))
    {
    printf("\n Give the distance in cm .........................: ");
    scanf("%lf",&xd);
    }
  
  double efield;
  if (phi>0)
    { // magentic deflection
      efield = 1.036427E-12*gamma*vel*vel*mass/z/rho;
    }
  else
    { // electrostatic deflection
      double gap;
      printf("\n Give the gap between the deflection plates in cm.: ");
      scanf("%lf",&gap);
      efield = 2.0*eion*log((rho+0.5*gap)/(rho-0.5*gap))/z;
    }
  double dwelltime = xd/vel/gamma;
  printf("\n The (motional) electric field strength is %5.2f kV/cm.",efield*0.001);
  printf("\n Do you want to change this value ? (y/n) ........: ");
  scanf("%s",answer);
  if ((!strcmp(answer,"y"))||(!strcmp(answer,"Y")))
    {
    printf("\n Give the (motional) E-field in V/cm .............: ");
    scanf("%lf",&efield);
    }
  printf("\n The dwell time in the field is %12.5g s\n",dwelltime);
  ncut = nf(efield,z);
  printf("\n Approximate cut-off quantum number nF =  %3d\n",ncut);
  printf("\n There are two models for the survival probabilities Ps,");
  printf("\n 1. a hard cut-off, i.e., Ps=1 for n<=ncut, Ps=0 for n>ncut");
  printf("\n 2. a soft cut-off calculated with FI rates ");
  printf("by Damburg and Kolosov\n");
  printf("\n Give nF (give 0 for soft cut-off) ..............:  ");
  scanf("%d",&ncut);
  softcut = ncut > 0 ? 0 : 1;

  strcpy(answer,"n");
  while (strcmp(answer,"y") && strcmp(answer,"Y"))
    {
    printf("\n Give nmax > ncut ...............................: ");
    scanf("%d",&nmax);

    if (nmax>nmaxfactorial/2.0) // nmaxfactorial is defined in factorial.h 
       {
       nmax = int(nmaxfactorial/2.0);
       printf("\n nmax exceeds nmaxfactoiral/2");
       printf(" which is defined in factorial.h");
       printf("\n The calculation of all necessary Clebsch-Gordan");
       printf("\n coefficients is not prossible.");
       printf("\n nmax set to nmaxfactorial/2 = %4d.\n",nmax);
       }
    double sigma_nmax=0.0, sigma_sum=0.0;
    for (int n=1;n<=nmax;n++)
       {
       sigma_nmax= 0.0;
       for (int l=0;l<n;l++)
         {
         sigma_nmax += sigmarrqm(1e-6,z,n,l);
         }
       sigma_sum += sigma_nmax;
       }
    printf("\n Near zero energy sigma(nmax) is %8.2g%% of",100*sigma_nmax/sigma_sum);
    printf(" the total RR cross section\n");
    printf("\n Is nmax o.k. ? (y/n) ............................: ");
    scanf("%s",&answer);
    } 


  // calculation of n,l-specific survival probabilities

  int n_one, n_zero;
  int psurvdim = nmax*(nmax+1)/2;
  double *psurv = new double[psurvdim];

  if (softcut)
    {
    calc_survival(dwelltime/gamma, efield, z, nmax, psurv, &n_one, &n_zero);
    }
  else
    {
    for (int n=1; n<=nmax; n++)
      {
      for (int l=0; l<n; l++) psurv[n*(n-1)/2+l] = n<=ncut ? 1 : 0;
      }
    n_one=ncut;
    n_zero=ncut+1;
    }

  // survival probabilities are 1 for n<=n_one and 0 for n>=n_zero
  // explicit calculations are only required for n_one < n < n_zero.

  printf("\n Now calculating hydrogenic transition rates and");
  printf(" decay probabilities  \n");

  RADRATE hydro(nmax,z,vel,x1/gamma,x2/gamma);
 
  double* sigma = new double[nmax];
  double* pd_list = new double[nmax];
  double* lt_list = new double[nmax];
  int* n_list = new int[nmax];
  int* l_list = new int[nmax];
  for(int i=0; i<nmax; i++) 
     {pd_list[i]=1.0; lt_list[i]=1.0; n_list[i]=1; l_list[i]=0;}

new_cascade:

  printf("\n\n Give number of cascade steps ( max if <0 ) ......: ");
  scanf("%d",&ncascstep);
  if (ncascstep<0) {ncascstep = nmax;}

  int read_previous = 0;
  printf("\n Multiply with values from a previous calculation?: ");
  scanf("%s",&answer);
  if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) )
    {
    read_previous = 1;

    printf("\n Give filename of previous calculation (*.fnl) ...: ");
    scanf("%s",&pfn);
    strcat(pfn,".fnl");

    char header[200];

    fin = fopen(pfn,"r");
    fgets(header,200,fin);
    printf("\n Header of previous calculation:\n %s \n",header);
 // overread next two lines
    fgets(header,200,fin);
    fgets(header,200,fin);
    }
      
  printf("\n Give filename for output (*.fnl, *.fm<i> i=1-3) .: ");
  scanf("%s",&filenameroot);

  strcpy(filename,filenameroot);
  strcat(filename,".fnl");
  fout = fopen(filename,"w");
  printf("\n Will write l-resolved results to %s.",filename);
  fprintf(fout,"z=%5.1f, ecool=%7.2f, x1=%7.2f cm, ",z,ecool,x1); 
  fprintf(fout,"x2=%7.2f cm, phi=%7.2f deg, ncascstep=%2d\n",
          x2,phi,ncascstep); 
  if (read_previous) {fprintf(fout," previous calculation: %s",pfn);}
  fprintf(fout,"\n   n    l            f  sigma(0 eV)      f*sigma\n");

  strcpy(filename,filenameroot);
  strcat(filename,".fn");
  fout2 = fopen(filename,"w");
  printf("\n Will write l-averaged survival probabilites to %s.",filename);

  strcpy(filename,filenameroot);
  strcat(filename,".fm1");
  fmatrix1 = fopen(filename,"w");
  printf("\n Will write l-resolved survival probabilities to %s.",filename);

  strcpy(filename,filenameroot);
  strcat(filename,".fm2");
  fmatrix2 = fopen(filename,"w");
  printf("\n Will write l-resolved survival probabilities weighted with RR cross section at 0 eV to %s.",filename);

  strcpy(filename,filenameroot);
  strcat(filename,".fm3");
  fmatrix3 = fopen(filename,"w");
  printf("\n Will write l-resolved decay probabilities to %s.",filename);

  printf("\n Now calculating detection probabilities \n");

  // calculate detection probabilities

  int n1, l1;
     
  for(n1=1; n1<=nmax; n1++)
    {
    double fn = 0.0;
    int n12 = (n1-1)*n1/2;
    // find maximum RR cross section per n and store cross sections
    // for latter use:

    double sigma_max = 0.0;
    for (l1=0; l1<n1; l1++)
      {
      sigma[l1] = sigmarrqm(1e-12,z,n1,l1);
      if (sigma[l1]>sigma_max) { sigma_max = sigma[l1];}
      }

    for(l1=0; l1<n1; l1++)
      {
      double tau = hydro.life(n1,l1);
      double pdecay = hydro.pdecay(n1,l1);
      double fnl = (1.0-pdecay)*psurv[n12+l1];
      
      if (((n1==1)||(n1==2))&&(l1==0))
	{
	  fnl = 1.0; // 1s and 2s states are assumed to be always detected
	}
      else
	{
	  pd_list[0] = pdecay; lt_list[0]  = tau; n_list[0] = n1; l_list[0] = l1;
      
	  int n2max;
	  if (ncascstep) {n2max = n1-1;}        // full calculation for cascades
	  else { n2max = (n1<=n_zero) ? n1-1 : n_zero;}
	  for(int n2=n2max; n2>0; n2--)
	    {
	      int n22 = (n2-1)*n2/2;
	      for(int l2=l1-1; l2<=l1+1; l2+=2)
		{
		  if ( ( l2<0) || (l2>=n2) ) continue;
		  double pcascdecay = 0.0, fcascade = 0.0;
		  if (ncascstep)
		    {
		      fcascade = cascade(0, ncascstep, n2, l2, n_zero, n_list, l_list,
                                pd_list, lt_list, hydro, psurv, &pcascdecay);
		    }
		  double br = hydro.branch(n1,l1,n2,l2);
		  fnl += br*((pdecay-pcascdecay)*psurv[n22+l2]+fcascade);
		} // end for (l2...)
	    } // end for (n2...)
	} // end else
      // read values from previous calculation

      int np,lp;
      double fnlp, sigmap, fsp;

      if (read_previous)
        {
        if (fscanf(fin,"%d %d %lf %lf %lf",&np,&lp,&fnlp,&sigmap,&fsp)==EOF) 
          {
          printf("\n !!! ATTENTION !!!");
          printf("\n Actual calculation not compatible with previous one!");
          printf("\n Program terminated at n = %4d.",n1-1);
          fclose(fin);
          fclose(fout);
          fclose(fmatrix1);
          fclose(fmatrix2);
          fclose(fmatrix3);
          break;
          }
        fnl*=fnlp;
        } // end if (read_previous)


      // output

      fprintf(fout,"%4d %4d %12.4G %12.4G %12.4G\n", 
                    n1,l1,fnl,sigma[l1],fnl*sigma[l1]/sigma_max);
      fprintf(fmatrix1,"%12.5g",fnl);
      fprintf(fmatrix2,"%12.5g",fnl*sigma[l1]/sigma_max);
      fprintf(fmatrix3,"%12.5g",pdecay);
      fn += (2.0*l1+1.0)*fnl;
      } // end for (l1...)

    fn /= n1*n1;
    fprintf(fout2,"%5d %12.4G\n",n1,fn);
  
    for (int l=n1; l<nmax; l++)
      {
      fprintf(fmatrix1,"%12.5g",0.0);
      fprintf(fmatrix2,"%12.5g",0.0);
      fprintf(fmatrix3,"%12.5g",0.0);
      }
    fprintf(fmatrix1,"\n");fprintf(fmatrix2,"\n");fprintf(fmatrix3,"\n");
     
    if ( (n1 % 10) == 0) { cout << "|"; } else { cout << "." ; } cout.flush();
                
    } // end for (n1...)
  if (read_previous) fclose(fin);
  fclose(fout);
  fclose(fout2);
  fclose(fmatrix1);
  fclose(fmatrix2);
  fclose(fmatrix3);
 
  printf("\n\n Another number of cascade steps? (y/n) : ");
  scanf("%s",&answer);
  if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) ) goto new_cascade;

  delete[] pd_list;
  delete[] lt_list;
  delete[] n_list;
  delete[] l_list;
  delete[] sigma;
  }

//////////////////////////////////////////////////////////////////////////

