/**
 * @file clebsch.cxx
 *
 * @brief Clebsch-Gordan coefficients, 3J, 6J, 9J symbols etc.
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: clebsch.cxx 375 2016-02-05 16:40:12Z iamp $
   @endverbatim
 *
 */
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include "clebsch.h"
#include "hydromath.h"

using namespace std;

inline int iround(double r) 
  { 
  return r < 0.0 ? (int) ceil(r-0.1) : (int) floor(r+0.1); 
  }

double delta(double a, double b, double c)
  { 
  return 0.5*(lnf(a+b-c)+lnf(a-b+c)+lnf(-a+b+c)-lnf(a+b+c+1));
  }



double SixJ(double j1, double j2, double j3, 
            double l1, double l2, double l3)
  {
  long double minus, sum1, sum2, s6j;
  int   k1, k2, k3, k4, k5, k6, k7, kmin, kmax, k;

  sum1 = delta(j1,j2,j3)+delta(j1,l2,l3)+delta(l1,j2,l3)+delta(l1,l2,j3);
  
  k1 = iround(j1+j2+j3);
  k2 = iround(j1+l2+l3);
  k3 = iround(l1+j2+l3);
  k4 = iround(l1+l2+j3);

  k5 = iround(j1+j2+l1+l2);
  k6 = iround(j2+j3+l2+l3);
  k7 = iround(j3+j1+l3+l1);

  kmin = k1;

  if (k2>kmin) kmin=k2;
  if (k3>kmin) kmin=k3;
  if (k4>kmin) kmin=k4;

  kmax = k5;
  if (k6<kmax) kmax = k6;
  if (k7<kmax) kmax = k7;
  
  s6j  =  0.0;
  for (k=kmin; k<=kmax; k++)
    {  
    minus = 1.0 - 2.0*abs(k % 2);
    sum2 = -lnf(k+1)+lnf(k5-k)+lnf(k6-k)+lnf(k7-k)+lnf(k-k1)
           +lnf(k-k2)+lnf(k-k3)+lnf(k-k4);
    s6j += minus*exp(sum1-sum2);
    }
  return double(s6j);
  } /* end SixJ */

double ThreeJold(double j1, double j2, double j3, 
              double m1, double m2, double m3)
  {
   long double sign, help, minus, sum1, sum2, s3j;
   int k1, k2, k3, k4, k5, kmin, kmax, k;

   if (iround(m1+m2+m3)!=0) return 0.0; 

   sum1 = delta(j1,j2,j3)+
           0.5*(lnf(j1-m1)+lnf(j1+m1)+lnf(j2-m2)
                +lnf(j2+m2)+lnf(j3-m3)+lnf(j3+m3));

   k1 = iround(j1+j2-j3);
   k2 = iround(j1-m1);
   k3 = iround(j2+m2);

   k4 = iround(j2-m1-j3);
   k5 = iround(j1+m2-j3);

   kmax = k1;
   if (k2 < kmax) kmax = k2;
   if (k3 < kmax) kmax = k3;

   kmin = 0;
   if (k4 > kmin) kmin = k4;
   if (k5 > kmin) kmin = k5;

   s3j  =  0.0;
   for (k=kmin; k<=kmax; k++)
     { 
     minus = 1.0 - 2.0*abs(k % 2);
     sum2  = lnf(k)+lnf(k1-k)+lnf(k2-k)+lnf(k3-k)+lnf(k-k4)+lnf(k-k5);
     s3j  += minus*exp(sum1-sum2);
     }
   sign = 1.0 - 2.0*abs(iround(j1-j2-m3) % 2);
   return double(sign*s3j);
  } /* end Three J*/

/////////////////////////////////////////////////////////////////////////
/**
 * @brief evaluation of the Three-J symbol in nested form 
 *
 * see I. Nasser and Y. Hahn, Phys. Rev. A 35, 2902 (1987)
*/
double ThreeJ(double j1, double j2, double j3, 
              double m1, double m2, double m3)
  {
   long double sign, help, minus, sum1, sum2, s3j;
   long double a, b, c, d, e, f;
   int k1, k2, k3, k4, k5, kmin, kmax, k;

   if (iround(m1+m2+m3)!=0) return 0.0; 

   sign = 1.0 - 2.0*abs(iround(j1-j2-m3) % 2);

   sum1 = delta(j1,j2,j3)+
           0.5*(lnf(j1-m1)+lnf(j1+m1)+lnf(j2-m2)
                +lnf(j2+m2)+lnf(j3-m3)+lnf(j3+m3));

   k1 = iround(j1+j2-j3);
   k2 = iround(j1-m1);
   k3 = iround(j2+m2);

   k4 = iround(j2-m1-j3);
   k5 = iround(j1+m2-j3);

   kmax = k1;
   if (k2 < kmax) kmax = k2;
   if (k3 < kmax) kmax = k3;

   kmin = 0;
   if (k4 > kmin) kmin = k4;
   if (k5 > kmin) kmin = k5;

   a = kmin;
   b = k1-kmin;
   c = k2-kmin;
   d = k3-kmin;
   e = kmin-k4;
   f = kmin-k5;

   minus = 1.0-2.0*abs(kmin % 2);
   long double *X = new long double[kmax-kmin+1];
   X[0] = lnf(a)+lnf(b)+lnf(c)+lnf(d)+lnf(e)+lnf(f);
   for (k=1;k<=kmax-kmin;k++)
       {
       X[k] = (b-k+1)*(c-k+1)*(d-k+1)/((a+k)*(e+k)*(f+k));
       }
   sum2  = 0.0;
   for (k=kmax-kmin;k>0;k--)
       {
       sum2 = (1.0-sum2)*X[k];
       }
   s3j = sign*(1.0-sum2)*minus*exp(sum1-X[0]);

   delete[] X;

   return double(s3j);
  } /* end Three J*/

double CG(double j1, double j2, 
          double m1, double m2, 
          double j, double m)
  {  
  int k;
  double minus;
  
  minus = 1.0 - 2.0*abs(iround(j1-j2+m) % 2);
  return minus*sqrt(2*j+1)*ThreeJ(j1,j2,j,m1,m2,-m);
  }

double Strength(double l1, double l2, double ml1, double dm)
  {
  double ml2, sl;

  ml2 = ml1+dm; 
  sl  = ThreeJ(l2,1.0,l1,-ml2,dm,ml1);
  return  (2*l1+1)*sl*sl;
  }


double NineJ( double ja1, double ja2, double ja3, 
              double jb1, double jb2, double jb3,
              double jc1, double jc2, double jc3)
{
    double ma1, ma2, ma3, mb1, mb2, mb3, mc1, mc2, mc3;

/*
    double sum = 0.0;
    for(ma1 = -ja1; ma1 <= ja1+0.1; ma1+=1)
    for(ma2 = -ja2; ma2 <= ja2+0.1; ma2+=1)
    for(ma3 = -ja3; ma3 <= ja3+0.1; ma3+=1)
    for(mb1 = -jb1; mb1 <= jb1+0.1; mb1+=1)
    for(mb2 = -jb2; mb2 <= jb2+0.1; mb2+=1)
    for(mb3 = -jb3; mb3 <= jb3+0.1; mb3+=1)
    for(mc1 = -jc1; mc1 <= jc1+0.1; mc1+=1)
    for(mc2 = -jc2; mc2 <= jc2+0.1; mc2+=1)
    for(mc3 = -jc3; mc3 <= jc3+0.1; mc3+=1)
    {
    printf("%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n",
	   ma1,ma2,ma3,mb1,mb2,mb3,mc1,mc2,mc3);
	sum +=
	    ThreeJ(ja1,ja2,ja3,ma1,ma2,ma3)*
	    ThreeJ(jb1,jb2,jb3,mb1,mb2,mb3)*
	    ThreeJ(jc1,jc2,jc3,mc1,mc2,mc3)*
	    ThreeJ(ja1,jb1,jc1,ma1,mb1,mc1)*
	    ThreeJ(ja2,jb2,jc2,ma2,mb2,mc2)*
	    ThreeJ(ja3,jb3,jc3,ma3,mb3,mc3);
    }

*/
    
    double jmax = jc2+jc3;
    if (jb1+jb3 > jmax) jmax = jb1+jb3;
    if (ja1+ja2 > jmax) jmax = ja1+ja2;
    double jmin = fabs(jc2-jc3);
    if (fabs(jb1-jb3) < jmin) jmin = jb1-jb3;
    if (fabs(ja1-ja2) < jmin) jmin = ja1-ja2;
 
    double sum = 0;
    for (double j = jmin; j<= jmax+0.1; j += 1.0)
    {
	double sign = 1-2*(int(2*j+0.1)%2);
	//printf(" %4.1f %4.1f\n",j,sign);
        sum += sign*(2*j+1)
	    * SixJ(ja1,jb1,jc1,jc2,jc3,j)
	    * SixJ(ja2,jb2,jc2,jb1,j,jb3)
	    * SixJ(ja3,jb3,jc3,j,ja1,ja2);
    }
    return sum;
}


void testCG(void)
  {
  float j1, j2, j, m1, m2, m, res;  

  printf("Give J1, M1 : ");scanf("%f %f",&j1,&m1);
  printf("Give J2, M2 : ");scanf("%f %f",&j2,&m2);
  printf("Give J , M  : ");scanf("%f %f",&j,&m);
  res = CG(j1,j2,m1,m2,j,m);
  printf("\n");
  printf("Clebsch Gordan coefficient  ");
  printf("<%4.1f,%4.1f,%4.1f,%4.1f| %4.1f,%4.1f> = %11.8f\n",
           j1,j2,m1,m2,j,m,res);
  printf("squared %16.8f\n\n",res*res);
  }

void StarkCG(void)
  {
  double res,n,k,m,m1,m2,j1,j2,j;
  int l;
  printf("Give n, k, m : ");
  scanf("%lf %lf %lf",&n,&k,&m);
  printf("\n");
  j1 = 0.5*(n-1.0);
  j2 = j1;
  m1 = 0.5*(m-k);
  m2 = 0.5*(m+k);
  for (l=int(m+0.1);l<n-0.1;l++)
    {
    j = l;
    res = CG(j1,j2,m1,m2,j,m);
    printf("l = %3d: CG = %11.8f, |CG|^2 = %11.8f\n",l,res,res*res);
    }
  printf("\n");
  }

void test3J(void)
  {
  float j1, j2, j3, m1, m2, m3, res;  

  printf("Give J1, J2, J3 : ");
  scanf("%f %f %f",&j1,&j2,&j3);
  printf("Give M1, M2, M3 : ");
  scanf("%f %f %f",&m1,&m2,&m3);
  res = ThreeJ(j1,j2,j3,m1,m2,m3);
  printf("\n");
  printf("           | %4.1f  %4.1f  %4.1f |\n",j1,j2,j3);
  printf("3-J Symbol |                  | = %11.8f\n",res);
  printf("           | %4.1f  %4.1f  %4.1f |\n",m1,m2,m3);
  printf("\n");
  }

void test6J(void)
  {
  float j1, j2, j3, l1, l2, l3, res;  

  printf("Give J1, J2, J3 : ");
  scanf("%f %f %f",&j1,&j2,&j3);
  printf("Give L1, L2, L3 : ");
  scanf("%f %f %f",&l1,&l2,&l3);
  res = SixJ(j1,j2,j3,l1,l2,l3);
  printf("\n");
  printf("           | %4.1f  %4.1f  %4.1f |\n",j1,j2,j3);
  printf("6-J Symbol |                  | = %11.8f\n",res);
  printf("           | %4.1f  %4.1f  %4.1f |\n",l1,l2,l3);
  printf("\n");
  }

void test9J(void)
  {
  float ja1, ja2, ja3;
  float jb1, jb2, jb3;
  float jc1, jc2, jc3;
  float res;

  printf("Give Ja1, Ja2, Ja3 : ");
  scanf("%f %f %f",&ja1,&ja2,&ja3);
  printf("Give Jb1, Jb2, Jb3 : ");
  scanf("%f %f %f",&jb1,&jb2,&jb3);
  printf("Give Jc1, Jc2, Jc3 : ");
  scanf("%f %f %f",&jc1,&jc2,&jc3);
  res = NineJ(ja1,ja2,ja3,jb1,jb2,jb3,jc1,jc2,jc3);
  printf("\n");
  printf("           | %4.1f  %4.1f  %4.1f |\n",ja1,ja2,ja3);
  printf("9-J Symbol | %4.1f  %4.1f  %4.1f | = %11.8f\n",jb1,jb2,jb3,res);
  printf("           | %4.1f  %4.1f  %4.1f |\n",jc1,jc2,jc3);
  printf("\n");
  }

void testStrength(void)
  {
  int l1, l2, m, dm;
  double sum, s;

  printf(" Give l->l' : ");
  scanf("%d %d",&l1,&l2);
  for (m = -l1; m <= l1; m++)
    {
    sum = 0.0;
    printf(" m=%2d",m);
    for (dm = -1; dm <= 1; dm++)
      {
      s = Strength(l1,l2,m,dm);
      sum += s;
      printf("%10.4f",s);
      }
    printf("%10.4f\n",sum);
    }
  }
