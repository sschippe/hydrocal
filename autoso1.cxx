/**
 * @file autoso1.cxx
 *
 * @brief Processing data from autostructure *.o1 output files
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: autoso1.cxx 374 2016-02-05 14:48:41Z iamp $
 @endverbatim
 *
 */

#include <stdio.h>
#include <math.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include "autoso1.h"
#include "lifetime.h"
#include "fele.h"

using namespace std;

AUTOSO1::~AUTOSO1()
{
    delete[] AA_CF1;
    delete[] AA_LV1;
    delete[] AA_W1;
    delete[] AA_CF2;
    delete[] AA_LV2;
    delete[] AA_rate;
    delete[] AA_energy;
    delete[] LV_K;
    delete[] LV_LV;
    delete[] LV_T;
    delete[] LV_S;
    delete[] LV_L;
    delete[] LV_J;
    delete[] LV_CF;
    delete[] LV_energy;
    delete[] AR_CF1;
    delete[] AR_LV1;
    delete[] AR_W1;
    delete[] AR_CF2;
    delete[] AR_LV2;
    delete[] AR_W2;
    delete[] AR_rate;
    delete[] AR_energy;
}

int AUTOSO1::open_file(char* filename)
{
    strcat(filename,".o1");
    fo1.open(filename,fstream::in);
    if (!fo1.is_open())
    {
	cout << "ERROR:: File " << filename << " does not exist!" << endl;
        return 0;
    }
    return 1;
}

void AUTOSO1::close_file(void)
{
    fo1.close();
}

int AUTOSO1::read_nl(int nRyd, int lRyd)
{
    const int AAdim=10000;
    const int LVdim=1000;
    const int ARdim=10000;

    AA_CF1 = new int[AAdim];
    AA_LV1 = new int[AAdim];
    AA_W1  = new int[AAdim];
    AA_CF2 = new int[AAdim];
    AA_LV2 = new int[AAdim];
    AA_rate = new double[AAdim];
    AA_energy = new double[AAdim];
    LV_K = new int[LVdim];
    LV_LV = new int[LVdim];
    LV_T = new int[LVdim];
    LV_S = new int[LVdim];
    LV_L = new int[LVdim];
    LV_J = new int[LVdim];
    LV_CF = new int[LVdim];
    LV_energy = new double[LVdim];
    AR_CF1 = new int[ARdim];
    AR_LV1 = new int[ARdim];
    AR_W1  = new int[ARdim];
    AR_CF2 = new int[ARdim];
    AR_LV2 = new int[ARdim];
    AR_W2 = new int[ARdim];
    AR_rate = new double[ARdim];
    AR_energy = new double[ARdim];
 
   
    char line[256], Cdummy;
    int nR, lR, nlevel, NucCharge, Nelectron;

    int n2=0, n3=0, n4=0, where=0;

//  read first line
    fo1.getline(line,256);
    while (!fo1.eof())
    {
	istringstream sline(line,istringstream::in);
        if (strncmp(line,"  NV=",5)==0)
	{ // check for next Rydberg electron
            if (where>1) 
	    {   //cout << " scan of file finished" << endl;
		break;
	    }
            where = 0;
            sline.seekg(5);
            sline >> nR;
            sline.seekg(15);
            sline >> lR;
            fo1.getline(line,256);
            if ( (nR==nRyd) && (lR==lRyd) )
	    {
		where = 1;
		cout << " Rydberg state nl: "<< nRyd << " " << lRyd << endl;
                nRydberg=nR; lRydberg=lR; // store Rydberg quantum numbers internally
	    }
 	}
        else if (strncmp(line,"        I-S            C-S",26)==0)
	{ // check for autoionization transition data
	    if (where>0) 
	    {
		where = 2;
                sline.seekg(66);
                sline >> NucCharge;
                sline.seekg(73);
                sline >> Nelectron;
		cout << "      nuclear charge: " << NucCharge << endl;
                cout << " number of electrons: " << Nelectron << endl;
                Zeff = NucCharge-Nelectron+1;
	    }
            fo1.getline(line,256);
	}
        else if (strncmp(line,"   NLEVEL=",10)==0)
	{ // check for level data
	    if (where>0) 
	    { 
		where = 3;
		sline.seekg(10);
		sline >> nlevel;
		cout << " number of levels: " << nlevel << endl;
	    }
            fo1.getline(line,256);
	}
        else if (strncmp(line,"        I-S            G-S",26)==0)
	{ // check for radiative transition data
	    if (where>0) 
	    {
		where = 4;
	    }
            fo1.getline(line,256);
	}
        else if (strcspn(line,"1234567890")<strlen(line))
	{ // if line containes any numbers
	    if (where==2)
	    { // read autoionization transition data
		sline >> AA_CF1[n2] >> AA_LV1[n2] >> AA_W1[n2];
		sline >> AA_CF2[n2] >> AA_LV2[n2] >> Cdummy;
		sline >> AA_rate[n2] >> AA_energy[n2];
		//cout << AA_CF1[n2] << " " << AA_LV1[n2] << " " << AA_W1[n2];
		//cout << " " << AA_CF2[n2] << " " << AA_LV2[n2];
		//cout << " " << AA_rate[n2] << " " << AA_energy[n2] << endl;
		n2++; 
                if (n2>AAdim) 
		{
		    cout << "ERROR: increase AAdim in autoso1.cxx" << endl;
		    return 0;
		} 
	    }
	    else if (where==3)
	    { // read level data
                sline >> LV_K[n3] >> LV_LV[n3] >> LV_T[n3] >> LV_S[n3] ;
                sline >> LV_L[n3] >> LV_J[n3] >> LV_CF[n3]; 
                sline >> LV_energy[n3] ;
                //cout << LV_K[n3] << " " << LV_LV[n3] << " " << LV_T[n3];
                //cout << " " << LV_S[n3] << " " << LV_L[n3] << " ";
                //cout << LV_J[n3] << " " << LV_CF[n3] << " ";
                //cout <<  LV_energy[n3] << endl;
		n3++;
                if (n3>LVdim)
		{
		    cout << "ERROR: increase LVdim in autoso1.cxx" << endl;
		    return 0;
		} 
	    }
	    else if (where==4)
	    { // read radiative transition data
		sline >> AR_CF1[n4] >> AR_LV1[n4] >> AR_W1[n4];
		sline >> AR_CF2[n4] >> AR_LV2[n4] >> AR_W2[n4];
		sline >> AR_rate[n4] >> AR_energy[n4];
		//cout << AR_CF1[n4] << " " << AR_LV1[n4] << " " << AR_W1[n4];
		//cout << " " << AR_CF2[n4] << " " << AR_LV2[n4] << " ";
                //cout << AR_W2[n4];
		//cout << " " << AR_rate[n4] << " " << AR_energy[n4] << endl;
		n4++;
                if (n4>ARdim) 
		{
		    cout << "ERROR: increase ARdim in autoso1.cxx" << endl;
		    return 0;
		} 
	    } 
	}
        // read next line
	fo1.getline(line,256);
    } // end while(!fo1.eof())
    nAA = n2-1;
    nLV = n3;
    nAR = n4;
    return 1;
}


int AUTOSO1::summed_rates(int nminHydro, int CF1)
{

    // sum of all hydrogenic radiative transition rate from n,l -> n'l' with n'>= nminHydro
    double  tau = hydrolife(Zeff,nRydberg,lRydberg,nminHydro);

    double *AAsum1 = new double[nLV];
    double *AAsum2 = new double[nLV];
    double *AAsum3 = new double[nLV];
    double *AAenergy1 = new double[nLV];
    double *AAenergy2 = new double[nLV];
    double *AAenergy3 = new double[nLV];
    double *ARsum = new double[nLV];
    int *AAcount1 = new int[nLV];
    int *AAcount2 = new int[nLV];
    int *AAcount3 = new int[nLV];
    int *ARcount = new int[nLV];
    int LV, i;
    double weight;
    for (i=0;i<nLV;i++)
    {
	AAsum1[i] = 0.0; 
	AAsum2[i] = 0.0;
	AAsum3[i] = 0.0;
        AAenergy1[i] = 0.0; 
        AAenergy2[i] = 0.0; 
        AAenergy3[i] = 0.0; 
	ARsum[i] = 1/tau;
        AAcount1[i]=0; 
        AAcount2[i]=0; 
        AAcount3[i]=0; 
        ARcount[i]=0; 
    }
    for (i=0;i<nAA;i++)
    {
	if (AA_CF1[i]==CF1)
	{ 
	    LV = AA_LV1[i];
            if (( (AA_CF2[i] ==-20) || (AA_CF2[i]==-17) ||
                  (AA_CF2[i] ==-14) || (AA_CF2[i]==-11) ||
                  (AA_CF2[i] ==-8) ) && (AA_energy[i] < 0.1) && (AA_rate[i]>0) )
	    { // rates for dielectronic capture from 3P0
                cout << LV << " " << AA_energy[i] << endl;
		AAsum1[LV-1] += fabs(AA_rate[i]);
                AAenergy1[LV-1] += AA_energy[i];
                AAcount1[LV-1]++;
	    }  
            else if (( (AA_CF2[i] ==-19) || (AA_CF2[i]==-16) ||
                  (AA_CF2[i] ==-13) || (AA_CF2[i]==-10) ||
                  (AA_CF2[i] ==-7) ) && (AA_rate[i]>0) )
	    { // rates for autoionization to the ground state
                //cout << AA_energy[i] << endl;
		AAsum2[LV-1] += AA_rate[i];
                AAenergy2[LV-1] += AA_energy[i];
                AAcount2[LV-1]++;
	    }
	    else
            {
		AAsum3[LV-1] += AA_rate[i];
                AAenergy3[LV-1] += AA_energy[i];
                AAcount3[LV-1]++;

	    }
	}
    }
    for (i=0;i<nAR;i++)
    {
	if ((AR_CF1[i]==CF1) && (AR_rate[i] > 0))
	{
	    LV = AA_LV1[i];
            ARsum[LV-1] += AR_rate[i];
	    ARcount[LV-1]++;
	}
    }
    cout << "   # 2J+1       Aa3P0   EDR3P0 ";
    cout << "       Aa1S0     E1S0         Arad";
    cout << "          EDR3P0      SDR      SRE  SRE/SDR  alphaDR" << endl ;
    for (i=0;i<nLV;i++)
    {
    if ( (AAsum1[i]>0.0))
        {
            for (int j=0; j<nLV; j++)
	    {
		if (LV_LV[j] == i) weight = LV_J[j]+1.0;
	    }
            cout.width(4);
            cout.precision(3);
	    cout << i << " ";
            cout.width(3);
            cout << weight << " ";
            cout.width(3);
            cout << AAcount1[i] << " ";
            cout.width(8);
            cout << AAsum1[i] << " ";
            cout.width(8);
            cout << 13.606*AAenergy1[i]/AAcount1[i] << " ";
            cout.width(3);
	    cout << AAcount2[i] << " ";
            cout.width(8);
            cout << AAsum2[i] << " ";
            cout.width(8);
            cout << 13.606*AAenergy2[i]/AAcount2[i] << " ";
//            cout.width(3);
//            cout << AAcount3[i] << " ";
//            cout.width(8);
//            cout << AAsum3[i] << " "; 
//            cout.width(8);
//            cout << 13.606*AAenergy3[i]/AAcount3[i] << " ";
            cout.width(3);
            cout << ARcount[i] << " ";
            cout.width(8);
	    cout << ARsum[i] << "        ";
            double vel = sqrt(2*AAenergy1[i]/511000)*3e10;
            double Edelta = 4*sqrt(0.69315*5E-5*AAenergy1[i]); // kTpar = 50 mueV
            double sum = AAsum1[i]+AAsum2[i]+ARsum[i];
	    double DRstrength = 4.95e-30*weight*AAsum1[i]*ARsum[i]/sum/2/AAenergy1[i];
	    double REstrength = 4.95e-30*weight*AAsum1[i]*AAsum2[i]/sum/2/AAenergy1[i];
            double DRalpha = DRstrength*vel*0.939/Edelta; // maximum of rate coefficient
            cout.width(8);
            cout << 13.606*AAenergy1[i]/AAcount1[i] << " ";
            cout.width(8);
            cout << DRstrength <<  " ";
            cout.width(8);
            cout << REstrength <<  " ";
            cout.width(8);
            cout << REstrength/DRstrength << " ";
            cout.width(8);
            cout << DRalpha << endl;
            
	} 
    }
    delete[] ARcount;
    delete[] AAcount3;
    delete[] AAcount2;
    delete[] AAcount1;
    delete[] ARsum; 
    delete[] AAenergy3;
    delete[] AAenergy2;
    delete[] AAenergy1;
    delete[] AAsum3;
    delete[] AAsum2;
    delete[] AAsum1;

	return 0;
}

void TestAutosO1(void)
{
    AUTOSO1 AO1;
    int nRyd, lRyd;

    char  filename[256];
    
    cout << "Give name of autostructure file (*.o1): ";
    cin >> filename;
    AO1.open_file(filename);

    cout << "Give Rydberg n and l: ";
    cin >> nRyd >> lRyd;
    AO1.read_nl(nRyd,lRyd);
    AO1.summed_rates(3,5);
    AO1.close_file();
}







