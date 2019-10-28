/**
 * @file kinema.cxx
 *
 * @brief Relativistic kinematics of merged electron-ion beams
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: kinema.cxx 446 2017-08-28 16:08:49Z iamp $
   @endverbatim
 *
 */

#include <math.h>
#include <stdio.h>
#include <string>
#include "kinema.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////

void kinema(void)
  {
  const double me = 0.511;
  const double amu = 931.494;
  FILE *fout;
  char fn[200];

  double E, Ee, Ei, A, Ec, Ed, Edmax, Edelta, theta;


  double mi,mu,mm,mm1,mm1sqr;
  double G, eps,ge,gi,sgige,bi,bisqr,gi1,ge1,gieps,sgieps;
  double c;
  double dEdG, dGdg, dGdeps, dEdg, dEdb, dEdeps;
  double Ex, dEdbx, dEdepsx, dbi, deps;

  printf("\n Give ion energy in MeV/u .............................: ");
  scanf("%lf",&Ei);
  printf("\n Give atomic mass in u ................................: ");
  scanf("%lf",&A);
  mi = A*amu;
  Ec = Ei/amu*me;
  gi = 1+Ei/amu;
  bisqr = 1.0-1.0/gi/gi;
  bi    = sqrt(bisqr);

  mu = mi*me/(mi+me);
  mm = me/mi;
  mm1 = 1.0+mm;
  mm1sqr =mm1*mm1;


  printf("\n Cooling energy:  %10.3f keV",1000*Ec);
  printf("\n ion gamma ....:  %10.6f",gi);
  printf("\n ion beta .....:  %10.6f",bi); 
  printf("\n\n Give max electron lab energy relative to cooling in eV: ");
  scanf("%lf",&Edmax);
  printf("\n Give step with in eV .................................: ");
  scanf("%lf",&Edelta);
  printf("\n Give relative uncertainty of ion beta ................: ");
  scanf("%lf",&dbi);
  printf("\n Give relative uncertainty of cathode voltage .........: ");
  scanf("%lf",&deps);

  printf("\n Give filename for output : ");
  scanf("%s",&fn);
  fout=fopen(fn,"w");
  
  theta = 0.0;
  c = cos(theta);
  printf("\n        Ed           E (eV)      ");
  printf("       DE(beta) (meV)          DE(Ed) (meV)"); 
  printf("\n                exact    approx.");
  printf("      exact    approx.      exact    approx.\n");
  for (Ed = 0; Ed <= Edmax; Ed+=Edelta)
    {
    Ee = Ed*1E-6+Ec;
    eps = Ed*1E-6/me;

    ge = 1.0+Ee/me;
    ge1 = ge*ge-1.0;
    gi1 = gi*gi-1.0;
    sgige = sqrt(gi1*ge1);
    gieps = (gi+0.5*eps)*eps/gi1;
    sgieps = sqrt(1.0+2.0*gieps);

    G = ge*gi-sgige*c;
    E = mi*mm1*(sqrt(1+2*mm*(G-1)/mm1sqr)-1.0);
    dEdG = me/mm1/sqrt(1+2*mm*(G-1)/mm1sqr);
    dGdeps = gi - ge/sgieps*c;
    dGdg = 2*gi+eps-2*gi*sgieps*c-(eps-2*gi*gieps)/sgieps*c;
    dEdb = fabs(dEdG*dGdg*gi*gi*gi*bisqr*dbi);
    dEdeps = fabs(dEdG*dGdeps*eps*deps);

//  approximate values
    Ex = 0.5*mu*(1-bisqr)/bisqr*eps*eps;
    dEdbx = 2*Ex/(1.0-bisqr)*dbi;
    dEdepsx = 2*Ex*deps;
    printf("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
           Ed,E*1E6,Ex*1E6,dEdb*1E9,dEdbx*1E9,dEdeps*1E9,dEdepsx*1E9);
    fprintf(fout,"%10.3f & %10.3f & %10.3f & %10.3f & %10.3f & %10.3f&  %10.3f \\ \n",
           Ed,E*1E6,Ex*1E6,dEdb*1E9,dEdbx*1E9,dEdeps*1E9,dEdepsx*1E9);
    }
  fclose(fout);
  }


///////////////////////////////////////////////////////////////////////////

double weizmass(int A, int Z, int Q)
  {
// returns nuclear mass for Q>=0 or 
// binding energy per nucleon for Q<0

  if ( (Z<1) || (A<1) || (Z>A) ) return 0.0;

  const double aV = 15.85e6;
  const double aS = 18.34e6;
  const double aC = 0.71e6;
  const double aA = 92.86e6;
  const double aP = 11.46e6;
  const double mH = 938.27231e6; // eV/c^2
  const double mn = 939.56563e6; // eV/c^2
  const double me = 0.51099906e6; // eV/c^2

  double A3,B,Tz;
  int N, sign, Zeven, Neven;

  N  = A - Z;
  Tz = 0.5*(Z - N);

  Zeven = Z % 2;
  Neven = N % 2;
  sign  = (Zeven==Neven) ? ( Zeven ? 1 : -1 ) : 0;

  A3 = exp(log(double(A))/3.0);
  B  = aV*A - aS*A3*A3 - aC*Z*Z/A3 - aA*Tz*Tz/A + sign*aP/sqrt(double(A));

  double wm = (Q>=0) ? Z*mH+N*mn-Q*me-B : B/double(A);
  
  return wm; 
  }

///////////////////////////////////////////////////////////////////////////

void TestWeizmass(void)
  {
  int a,z,q,mode,amax,amin;
  double wm;
  FILE *fout;
  char fn[200];

  const double malpha = 3727.409e6; // eV/c^2 
  printf("\n Give maximum mass number : ");
  scanf("%d",&amax);
  printf("\n Mode 0 : max(B(Z,A)/A) as function of A");
  printf("\n Mode 1 : B(Z,A)/A");
  printf("\n Mode 2 : m(Z,A)");
  printf("\n Mode 3 : m(Z,A)-m(Z-2,A-4)-m(alpha)");
  printf("\n Give mode ...............: ");
  scanf("%d",&mode);
  printf("\n Give filename for output : ");
  scanf("%s",&fn);

  fout=fopen(fn,"w");
  amin = (mode<3) ? 1 : 5;
  q = (mode>1) ? 0 : -1;
  for(a=amin; a<=amax; a++)
    {
    int zmax=1;
    double wmax = 0.0;
    for(z=amin; z<=amax; z++)
      {
      wm = weizmass(a,z,q);
      if (wm>wmax) { wmax=wm, zmax=z; }
      if (mode>2) 
         {
         double wm2 = weizmass(a-4,z-2,q);
         if (wm2>0) { wm -= wm2+malpha;} else {wm = -malpha;}
         }
      if (mode>0) fprintf(fout,"%12.6g ",wm/1e6);
      }
    if (mode>0) {fprintf(fout,"\n");}
    else 
      {
      double wzmax = weizmass(a,zmax,0);
      double ealpha = wzmax-weizmass(a-4,zmax-2,0)-malpha;
      double ecoul = -1.11e6*z*z*exp(-log(0.5*a)/3)/8.0;
      double efiss = weizmass(a/2,zmax/2,0);
      efiss = (efiss>0) ? wzmax-2*efiss : 0.0;
      fprintf(fout,"%5d %5d %12.6g %12.6g %12.6g %12.6g %12.6g\n",
              a,zmax,wmax/1e6,wzmax/1e6,ealpha/1e6,efiss/1e6,ecoul/1e6);
      }
    }
  fclose(fout);
  printf("\n Mass output  written to file %s.\n\n",fn);
  }

///////////////////////////////////////////////////////////////////////////
