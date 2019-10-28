/**
 * @file stripping.cxx
 *
 * @brief Formulae for equilibrium charge states resulting from stripping
 *
 * @author Stefan.Schippers@physik.uni-giessen.de
 *
 * @copyright GNU General Public License.
 *
 * $Id: stripping.cxx 375 2016-02-05 16:40:12Z iamp $
 *                                                        						
 * Modification Log:
 *     - v1.0, 2013-10-30 by Stefan Schippers
 */ 

#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;
////////////////////////////////////////////////////////////////////////////////////
/**
* @brief Mean charge after stripping in gas
*
* G. Schiwietz and J. P. Grande, NIMB 175-177 (2001) 125-131, Eqs. 1 & 2
*
* @param vp projectile velocity in atomic units
* @param Zp nuclear charge of projectile
* @param Zt nulcear charge of target
* @param x on exit, reduced scaling variable
*
* @return mean charge state
*/
double MeanChargeGas(double vp, double Zp, double Zt, double &x)
{
	double vZp = vp*pow(Zp,-0.52);
	x = pow(vZp*pow(Zt,0.03-0.017*vZp),1.0+0.4/Zp);
	double x3 = x*x*x;
	double x6 = x3*x3;
	return Zp*(376.0*x+x6)/(1428.0-1206.0*sqrt(x)+690.0*x+x6);
}

////////////////////////////////////////////////////////////////////////////////////
/**
* @brief Mean charge after stripping in solid
*
* G. Schiwietz and J. P. Grande, NIMB 175-177 (2001) 125-131, Eqs. 3 & 4
*
* @param vp projectile velocity in atomic units
* @param Zp nuclear charge of projectile
* @param Zt nulcear charge of target
* @param x on exit, reduced scaling variable
*
* @return mean charge state
*/
double MeanChargeSolid(double vp, double Zp, double Zt, double &x)
{
	double vZp = vp*pow(Zp,-0.52);
	x = pow(vZp*pow(Zt,-0.019*vZp)/1.68,1.0+1.8/Zp);
	double x2 = x*x;
	double x4 = x2*x2;
	return Zp*(12.0*x+x4)/(0.07/x+6.0+0.3*sqrt(x)+10.37*x+x4);
}


////////////////////////////////////////////////////////////////////////////////////
/**
* @brief Width of charge state distribution after stripping in solid
*
* G. Schiwietz and J. P. Grande, NIMB 175-177 (2001) 125-131, Eqs. 5 & 6
*
* @param Zp nuclear charge of projectile
* @param Zt nulcear charge of target
* @param qmean mean charge state
*
* @return mean charge state
*/
double MeanChargeWidth(double Zp, double Zt, double qmean)
{
	const double w = 0.7;
	double a = 0.37*pow(Zp,0.6);
	double f1 = sqrt((qmean+a)/qmean);
	double f2 = sqrt((Zp-qmean+a)/(Zp-qmean));
	return w/(pow(Zp,-0.27)*pow(Zt,0.035-0.0009*Zp)*f1*f2);
}


////////////////////////////////////////////////////////////////////////////////////
/**
* @brief Gaussian charge state distribution after stripping in solid
*
* @param q charge state
* @param qmean mean charge state
* @param Gaussian width
*
* @return charge state fraction
*/
double MeanChargeGauss(double q, double qmean, double width)
{
	const double fac = 0.39894228040143;  // 1/sqrt(2 pi)
	double arg = (q-qmean)/width;
	return fac/width*exp(-0.5*arg*arg);
}

////////////////////////////////////////////////////////////////////////////////////
/**
* @brief Calculates charge state distributions after stripping in solids and gases
*
* The user is prompted for input. Results are written to a worksheet.
*/
void stripping(void)
{
    const double muc2 = 931.494028;     // [MeV] energy equivalent of atomic mass unit
    const double clight = 137.035999074;  // inverse fine structure constant

	double Ep, Mp, Zp, Zt;
	int striptype;
	char answer;

	cout << endl;
	cout << " Calculation of Gaussian equilibrium charge state distributions" << endl;
	cout << " after stripping, using semi-empirical formulae of G. Schiwietz" << endl;
	cout << " and J. P. Grande [NIMB 175-177 (2001) 125-131]." << endl;
	cout << endl << " Give type of stripper, foil or gas (f/g): ";
    cin >> answer;
	if ((answer=='f') || (answer=='F'))
	{
		striptype = 1;
		cout << endl << " Give nuclear charge of foil target .......: ";
	}
	else
	{
		striptype = 0;
		cout << endl << " Give nuclear charge of gas target ........: ";
	}
	cin >> Zt;

	cout << endl << " Give nuclear charge of projectile ........: ";
	cin >> Zp;

	cout << endl << " Give atomic mass number of projectile ....: ";
	cin >> Mp;

	for(;;)
	{
		cout << endl << " Give projectile energy in MeV (0 quits) ..: ";
		cin >> Ep;
		if (Ep<0.1) break;
    
		double gamma = 1.0+Ep/(Mp*muc2);
		double vp = clight*sqrt(1.0-1.0/(gamma*gamma));

		double qmean, x;
		if (striptype==0)
		{
			qmean = MeanChargeGas(vp,Zp,Zt,x);
		}
		else
		{
			qmean = MeanChargeSolid(vp,Zp,Zt,x);
		}
		double width = MeanChargeWidth(Zp,Zt,qmean);

		cout << endl << " mean charge state : " << qmean;
		cout << endl << " distribution width: " << width << endl;
		int count = 0;
		double thres = 0.1;
		for (int q=0; q<=Zp+0.1; q++)
		{
			double frac = 100.0*MeanChargeGauss(q,qmean,width);
			frac = 0.1*floor(10*frac+0.5);
			if (frac>thres) {
				cout << " fq(" << setw(2) << q  << ") = ";
				cout << setw(5) << setprecision(1) << fixed  << frac << " %" << endl;
			}
		}
	} // end for
}
