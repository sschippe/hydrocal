/**
 * @file burgess.cxx
 *
 * @brief Burgess formula for atomic photoionization cross sections (untested)
 *
 * see A. Burgess and J. M. Seaton, MNRAS 120 (1960) 121--151
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: burgess.cxx 375 2016-02-05 16:40:12Z iamp $
   @endverbatim
 *
 */
#include <stdio.h>
#include <math.h>

using namespace std;


double g(double eps, double nu, double mu, int l1, int l2, double *cosphi)
  {
// Table III

  double  a[6] = {-0.147,-0.216,-0.120, -0.247,-0.117,-0.362};
  double  b[6] = {0.2515,-0.171, 0.600, -0.272, 1.170, 0.599};
  double  c[6] = {-0.078, 0.000, 0.000,  0.000, 0.000,-2.432};
  double aa[6] = { 0.310, 0.000, 0.362, -0.010, 0.321,-0.390};
  double bb[6] = { 0.000, 0.000, 0.0535,-0.019, 0.106, 0.050};
 
// Table IV

  double  G[6][12] = {{2.723, 2.095, 1.856, 1.718, 1.623, 1.553, 
                       1.498, 1.452, 1.414, 1.381, 1.352, 1.327},
                      {0.000, 1.028, 1.117, 1.152, 1.168, 1.175,
                       1.177, 1.176, 1.173, 1.170, 1.165, 1.161},
                      {0.000, 2.840, 2.264, 2.010, 1.856, 1.749,
                       1.666, 1.601, 1.546, 1.501, 1.461, 1.427},
                      {0.000, 0.000, 0.669, 0.818, 0.899, 0.952,
                       0.988, 1.014, 1.033, 1.047, 1.058, 1.065},
                      {0.000, 0.000, 3.000, 2.143, 2.139, 1.971,
                       1.854, 1.765, 1.694, 1.635, 1.585, 1.543},
                      {0.000, 0.000, 0.000, 0.468, 0.599, 0.704,
                       0.793, 0.868, 0.933, 0.991, 1.041, 1.085}};
                            
// Table V

  double gg[6][12] = {{1.754, 1.605, 1.591, 1.590, 1.591, 1.594, 
                       1.596, 1.599, 1.601, 1.603, 1.605, 1.607},
                      {0.000, 1.667, 1.667, 1.667, 1.667, 1.667,
                       1.667, 1.667, 1.667, 1.667, 1.667, 1.667},
                      {0.000, 1.574, 1.582, 1.589, 1.582, 1.587,
                       1.593, 1.598, 1.603, 1.608, 1.614, 1.618},
                      {0.000, 0.000, 1.819, 1.771, 1.741, 1.722,
                       1.707, 1.697, 1.688, 1.682, 1.676, 1.672},
                      {0.000, 0.000, 1.447, 1.535, 1.544, 1.549, 
                       1.556, 1.564, 1.573, 1.581, 1.589, 1.596},
                      {0.000, 0.000, 0.000, 1.850, 1.908, 1.918,
                       1.920, 1.921, 1.922, 1.924, 1.926, 1.928}};


// Table VI

  double  G01[8] = { 3.259, 2.976, 2.739, 2.527, 2.360, 2.244, 2.162, 2.095};
  double gg01[8] = { 1.850, 1.770, 1.701, 1.655, 1.632, 1.620, 1.612, 1.604};
  double xx01[8] = { 0.143, 0.085, 0.043, 0.011,-0.008,-0.020,-0.031,-0.041};

// Table VII

  double Gn10[11] = { 1.880, 1.500, 1.310, 1.180, 1.076, 0.999, 0.939,
                      0.894, 0.854, 0.820, 0.790};
  double  G10[11] = { 0.000, 0.670, 0.826, 0.911, 0.962, 0.999, 1.029,
                      1.058, 1.080, 1.100, 1.117};
  double gg10[11] = { 1.333, 1.515, 1.585, 1.630, 1.655, 1.667, 1.667,
                      1.667, 1.667, 1.667, 1.667};
  double xx10[11] = {-0.330,-0.321,-0.313,-0.306,-0.300,-0.295,-0.290,
                     -0.286,-0.281,-0.277,-0.273};

// Table VIII

  double Gn12[11] = { 5.690, 5.570, 5.020, 4.250, 3.438, 2.757, 1.095,
                      0.879, 1.886, 1.342, 0.735};
  double  G12[11] = { 0.000, 2.489, 3.174, 3.291, 3.075, 2.757, 2.512,
                      2.415, 2.386, 2.340, 2.251};
  double gg12[11] = { 2.340, 1.911, 1.703, 1.625, 1.624, 1.658, 1.675,
                      1.635, 1.593, 1.576, 1.597};
  double xx12[11] = { 0.650, 0.511, 0.389, 0.287, 0.210, 0.164, 0.1425,
                      0.131, 0.115, 0.0945, 0.073};
  double bb12[11] = { 0.079, 0.069, 0.054, 0.038, 0.029, 0.035, 0.053,
                      0.068, 0.068, 0.060, 0.050};


  const double pi = 3.1415926536;
  double en, en1, enn, enn1, nn, nn5, xx, ggn, bbn, Gn, Gns, cp, sign;
  int n, n5, i;
  
  i = 2*l1 + (l2-l1+1)/2 - 1;
  
  n   = int(nu); nn = nu-n; 
  en  = eps*nu; en1 = 1.0+en; 
  enn = en*nu; enn1 = 1.0+enn;

  if ( nu >= (l1+2) )
    {
    xx = a[i]+(b[i]+c[i]/nu)/nu + en/en1*aa[i] + enn/enn1*bb[i];
    ggn = gg[n-1][i]*(1-nn) + gg[n][i]*nn; 
    Gn  =  G[n-1][i]*(1-nn) +  G[n][i]*nn; 
    Gns = Gn/sqrt(1.0);
    }
  else
    {
    nn5 = nu*5; // interpolation step size 0.2
    n5  = int(nn5); nn5 -= n5;
    switch (i)
      {
      case 0: xx  = xx01[n5-3]*(1-nn5) + xx01[n5-2]*nn5;
	      xx += 0.310*en/en1;
              ggn = gg01[n5-3]*(1-nn5) + gg01[n5-2]*nn5;
               Gn =  G01[n5-3]*(1-nn5) +  G01[n5-2]*nn5;
              Gns = Gn/sqrt(nu-1);
              break;

      case 1: xx  = xx10[n5-5]*(1-nn5) + xx10[n5-4]*nn5;
              ggn = gg10[n5-5]*(1-nn5) + gg10[n5-4]*nn5;
              Gns = Gn10[n5-5]*(1-nn5) + Gn10[n5-4]*nn5;
              break;
  
      case 2: xx  = xx12[n5-5]*(1-nn5) + xx12[n5-4]*nn5;
              bbn = bb12[n5-5]*(1-nn5) + bb12[n5-4]*nn5;
              xx += 0.362*en/en1 + enn/enn1*bbn;
              ggn = gg12[n5-5]*(1-nn5) + gg12[n5-4]*nn5;
              Gns = Gn12[n5-5]*(1-nn5) + Gn12[n5-4]*nn5;
              break;

      default: printf("error in g!\n\n"); return 0.0; break;

      } // end switch
    } // end else 
  //printf("%6.3f %6.3f %6.3f %6.3f\n",nu,Gns,ggn,xx);

  sign = 2*(l1%2-0.5);
  cp = cos(pi*(nu+mu+xx)); *cosphi = cp;
  return sign*Gns*exp(-ggn*log(enn1))*cp;
  }

double sigmaPI(double e, double z, double nu, double mu, int l, int w,
               double *glm, double *glp, double *cosphim, double *cosphip)
  {
  double gm, gp, Cm, Cp, phip, phim, z2, Inl;

  z2 = z*z; Inl = z2/nu/nu;

  *cosphim = 0.0;
  gm = (l==0) ? 0 : g(e, nu, mu, l, l-1, cosphim);
  gp = g(e, nu, mu, l, l+1, cosphip);

  *glm = gm;
  *glp = gp;

  switch (l)
    {
    case 0: Cm = 0.0; Cp = w; break;   // s^w
    
    case 1: Cm = w/3.0; Cp = 2*w/3.0; break; // p^w, w<=2
         
    default: Cm= 0.0; Cp = w; break;

    } // end switch(l)


  return 8.559e-19*(Inl+e*z2)/(Inl*Inl)*(Cm*gm*gm + Cp*gp*gp);
  
  }


void BurgessPIxsec(void)
  {
  double e, eI, z, nu, mu, mu0, mu1, emax, sigma, gm, gp, cm, cp;
  int n, l, w;
  char filename[200];
  FILE *fout;

  printf(" Outer shell photoionization of ions with one open nl^w shell\n");
  printf("\n Give main quantum number n ...............: ");
  scanf("%d",&n);
  printf("\n Give l (0,1) .............................: ");
  scanf("%d",&l);
  printf("\n Give w (1,2) .............................: ");
  scanf("%d",&w);
  printf("\n Give ionization energy in eV .............: ");
  scanf("%lf",&eI);
  printf("\n Give charge state after ionization .......: ");
  scanf("%lf",&z);
  printf("\n Give quantum defect of continuum electron : ");
  scanf("%lf",&mu0); 
  if (mu0!=0.0)
    {
    printf("\n The quantum defect is modelled as mu0+E*mu1");
    printf("\n Give coefficient of linear term ..........: ");
    scanf("%lf",&mu1);
    }
  printf("\n Give max electron energy in Rydberg ......: ");
  scanf("%lf",&emax);
  printf("\n Give filename for output .................: ");
  scanf("%s",&filename);

  nu = z*sqrt(13.606/eI);
  printf("\n The effective quantum number is ........ nu = %7.4f\n",nu);

  fout = fopen(filename,"w");printf("\n");
  for (e=0.0; e<=1.02*emax; e+=0.05*emax)
    {
    mu = mu0+e*mu1;
    sigma = sigmaPI(e,z,nu,mu,l,w,&gm,&gp,&cm,&cp);
    fprintf(fout,"%12.4g %12.4g %12.4g %6.3f %12.4g %6.3f\n",
           e,sigma,gm,cm,gp,cp);
    printf("%12.4g %12.4g %12.4g %6.3f %12.4g %6.3f\n",
           e,sigma,gm,cm,gp,cp);
    }
  fclose(fout);
  }











