/**
 * @file BetheBloch.cxx
 *
 * @brief Stopping power calculations with the Bethe-Bloch formula
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: BetheBloch.cxx 375 2016-02-05 16:40:12Z iamp $
   @endverbatim
 *
 */
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>

using namespace std;
const int maxcomposits=10;
const double atomic_mass=931.49 ; //MeV
const double electron_mass=0.511;// MeV

/**
 * @brief Bethe-Bloch formula for the stopping of a heavy charged particle in dense matter
 *
 * @param ep projectile energy in MeV
 * @param zp projectile nuclear charge
 * @param ap projectile atomic mass number
 * @param zt target nuclear charge
 * @param at target atomic mass number
 * @param it target first ionization potential in eV
 *
 * @return stopping power in Me cm²/g
 */

double Stopping(double ep, double  zp, double ap, double zt, double at, double it){
    double gamma= 1+(ep/ap)/atomic_mass;
    double gammasqr=gamma*gamma;
    double betasqr= 1.0 - 1.0/gammasqr;
    double term= 2000000*electron_mass*gammasqr*betasqr/it;
    if(term<1){return 0;} else{ return 0.307*(zt/at)*zp*(zp/betasqr)*(log(term)-betasqr);}
};

/**
 * @brief Bragg's rule for the stopping power of a compound
 *
 * @param ep projectile energy in MeV
 * @param zp projectile nuclear charge
 * @param ap projectile atomic mass number
 * @param zt target nuclear charge (vector with one entry for each atomic contituent)
 * @param at target atomic mass number (vector with one entry for each atomic contituent)
 * @param it target first ionization potential in eV (vector with one entry for each atomic contituent)
 * @param wt weight of each constituent (vector with one entry for each atomic contituent)
 * @param stp on exit stopping power per constituent (vector with one entry for each atomic contituent)
 * @param number_of_composits number of atomic constituents
 *
 * @return stopping power of the compound in MeV cm² / g 
 */
double Compound(double ep, double zp, double ap,vector <double> zt,vector <double> at,vector <double> it, vector <double> wt, vector <double> &stp, int number_of_composits){
    double sum=0;
    for (int n=1;n<=number_of_composits;n++){
	stp[n] = Stopping(ep,zp,ap,zt[n],at[n],it[n]);
        sum+= wt[n]*stp[n];
    }
    return sum;
};


/**
 * @brief User input of target properties
 *
 * @param zt on exit target nuclear charge (vector with one entry for each atomic contituent)
 * @param at on exit target atomic mass number (vector with one entry for each atomic contituent)
 * @param ct on exit number of atoms per constituent (vector with one entry for each atomic contituent)
 * @param it on exit target first ionization potential in eV (vector with one entry for each atomic contituent)
 * @param wt weight of each constituent (vector with one entry for each atomic contituent)
 * @param density on exit density of the target
 *
 * @return number of constituents in compound 
 */
int TargetInput(vector <double> &zt, vector <double> &at, vector <double> &ct, vector <double> &it, vector <double> &wt, double &density){
    int choice=0,n=0,nt=0;
    double sum=0;
    bool weiter=false;
    do {
        cout<<endl;
        cout<<" Target selection"<<endl;
        cout<<"   1: manual input"<<endl;
        cout<<"   2: water"<<endl;
        cout<<"   3: mylar"<<endl;
        cout<<"   4: havar"<<endl;
        cout<<"   5: aluminium"<<endl<<endl;
        cout<<" Make your choice (0 quits )               : ";
        cin>>choice;
        switch (choice) {
          case 0 :    return 0;
          case 1 :    cout<<endl<<endl<<" Manual input of target components (maximal 10)"<<endl;
                      for (nt=1;nt<=maxcomposits;nt++) {
                          cout<<    "   Give target component No "<<nt<<endl;
                          cout<<    "    nuclear charge (0 quits) .......: ";
                          cin>>zt[nt];
                          cout<<endl;
                          if (zt[nt]>0.5){
                              cout<<"   Give target atomic masse  .......: ";
                              cin>>at[nt];
                              cout<<endl;
                              cout<<"   Give target atom abundance ......: ";
                              cin>>ct[nt];
                              cout<<endl;
                          } else {nt--; break;}
                      }
                      if (!nt) break;
                      cout<<endl<<  "   Give target density in g/cm^3 ...: ";
                      cin>>density;
                      cout<<endl;
                      break;
          case 2 :    nt=2;
                      density=1;
                      zt[1]=1;
                      at[1]=1.008;
                      ct[1]=2;
                      zt[2]=8;
                      at[2]=16;
                      ct[2]=1;
                      break;
          case 3 :    nt=3;
                      density=1.395;
                      zt[1]=6;
                      at[1]=12.01;
                      ct[1]=10;
                      zt[2]=1;
                      at[2]=1.008;
                      ct[2]=8;
                      zt[3]=8;
                      at[3]=16;
                      ct[3]=4;
                      break;
          case 4 :    nt=6;
                      density=8.3;
                      zt[1]=27;
                      at[1]=58.9332;
                      ct[1]=42.5;
                      zt[2]=24;
                      at[2]=51.9961;
                      ct[2]=20;
                      zt[3]=25;
                      at[3]=54.93805;
                      ct[3]=1.6;
                      zt[4]=42;
                      at[4]=95.94;
                      ct[4]=2.0;
                      zt[5]=28;
                      at[5]=58.69;
                      ct[5]=13;
                      zt[6]=74;
                      at[6]=183.85;
                      ct[6]=2.8;
                      break;
          case 5 :    nt=1;
                      density=2.7;
                      zt[1]=13;
                      at[1]=26.98;
                      ct[1]=1;
                      break;
           default :  nt=0;
        }
        choice=0;
    } while (!nt);
    for (n=1;n<=nt;n++) { sum+=at[n]*ct[n]; }
    cout<<" Target parameters are:"<<endl;
    for (n=1;n<=nt;n++) {
        wt[n]=at[n]*ct[n]/sum;
        if (zt[n]<13) {
            it[n]=12*zt[n]+7;
        } else {
            it[n]=9.79*zt[n]+58.8*exp(-0.19*log(zt[n]));
        }
        printf("   z%i: %7.3f  a%i: %7.3f  c%i: %7.3f  w%i: %7.3f  i%i: %7.3f\n",
                     n, zt[n],   n, at[n],   n, ct[n],   n, wt[n],   n, it[n] );
    }
    return nt;
}
/**
 * @brief Entry to Bethe-Bloch stopping power calculation
 */
int BetheBloch()
{
    vector<double> zt(maxcomposits,0);
    vector<double> at(maxcomposits,0);
    vector<double> ct(maxcomposits,0);
    vector<double> wt(maxcomposits,0);
    vector<double> it(maxcomposits,0);
    vector<double> stp(maxcomposits,0);
    int number_of_composits=0;
    int n=0;
    double zp=0, ap=0, density=0, sum=0, ep=0, ep0, stop=0, x=0, dx=0;
    char filename[200];
    number_of_composits=TargetInput(zt,at,ct,it,wt,density);
    if (number_of_composits==0) return 0;
    cout<<endl;
    cout<<endl<<"   Give nuclear charge of projectile .: ";
    cin>>zp;
    cout<<endl<<"   Give atomic mass of projectile in u: ";
    cin>>ap;
    cout<<endl<<"   Give initial energy in MeV ........: ";
    cin>>ep;
    ep0 = ep;
    cout<<endl;
    cout<<endl<<"   Give step width in micrometers"<<endl;
    cout<<      "    (Don't use to small step width!)..: ";
    cin>>dx;
    dx=dx*1e-4;
    cout<<endl;
    cout<<endl<<"   Give filename for output ..........: ";
    scanf("%s",&filename);
    cout<<endl<<endl;
    FILE* fout=fopen(filename,"w");
    fprintf(fout, "       x            E         DELTA E        -dE/dx     ");
    for (int n=1; (n<=number_of_composits) && (number_of_composits>1); n++)
    {
	fprintf(fout,"      q%1d=%5.0f ",n,zt[n]);
    }
    fprintf(fout,"\n");
    fprintf(fout, "        [mu]        [MeV]        [keV]     [MeV cm^2/g] ");
    for (int n=1; (n<=number_of_composits) && (number_of_composits>1); n++)
    {
	fprintf(fout,"      w%1d=%5.3f ",n,wt[n]);
    }	
    fprintf(fout,"\n");
    stop=1;
    while((stop>0) && (ep>0.05))
    {
	double eloss = (ep0-ep)*1.0E3; // energy loss
        stop=Compound(ep,zp,ap,zt,at,it,wt,stp,number_of_composits);
        fprintf(fout,"%12.6f %12.6f %14.6E %14.6E",x*1e4,ep,eloss,stop);
        for (int n=1; (n<=number_of_composits) && (number_of_composits>1); n++)
	{
		fprintf(fout," %14.6E",stp[n]);
	}
	fprintf(fout,"\n");
        ep -= stop*density*dx;
        x+=dx;
    }
    fclose(fout);
    return 0;
}

