/**
 * @file peakform.cxx
 *
 * @brief LaTeX formatting of Origin peak fitting results (obsolete?)
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: peakform.cxx 375 2016-02-05 16:40:12Z iamp $
   @endverbatim
 *
 */
#include <strstream>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

using namespace std;

void peakformat(void)
{
	char fn[200],fn_in[200],fn_out[200], line[200];
	cout << " Give name of file containing FanoVoigt fit results (*.dat): ";
	cin >> fn;
    strcpy(fn_in,fn); strcat(fn_in,".dat");
	strcpy(fn_out,fn); strcat(fn_out,".tex");

	cout << '\n' << " reading file " << fn_in << "\n\n";
	cout.flush();
    
	ifstream fin(fn_in);
	ofstream fout(fn_out);
    istrstream buf(line,200);
	char symb[20], teststr[3];
	double val,err,fval,ferr;
	int startflag=0, count=0, i;

	fout << "\\begin{tabular}{ccccc}\n";
    fout << " \\hline \\hline \n";
	fout << " resonance & $E_\\mathrm{res}$ & $Q$ & $\\Gamma$ & $A$ \\\\\n";
    fout << " \\hline \n";
	while (!fin.eof())
	{
		buf.clear();
		fin.getline(line,200);
	    for (int i=0;i<2;i++) teststr[i]=line[i]; teststr[2] = '\0';
		if ((!startflag) && (!strcmp(teststr,"e1"))) startflag=1;
		if (startflag)
		{
            if (!strcmp(teststr,"--")) break;
			if (count<4)
			{
				fout << "  &  ";
			}
			else
			{
				fout << " \\\\ \n  &  ";
				count = 0;
			}
			count++;
            buf.seekg(0);
			buf >> symb >> val >> err;
			cout.width(10);
			cout << symb;
	    	cout << " = ";
			cout.width(12);
			cout << val <<" +- ";
			cout.width(12);
			cout << err;
			cout <<" formatted: ";
            if (err==0.0)
			{
				cout << val <<"(-)\n";
				fout << val <<"(-)";
				continue;
			}
			int sign = val < 0 ? -1 : 1;
            int lval = val != 0 ? int(log10(fabs(val))) : 0;
			int lerr = err != 0 ? int(log10(fabs(err))) : 0;
		    lval = lval > 0 ? lval+1 : lval-1;
		    lerr = lerr > 0 ? lerr+1 : lerr-1;
		    if (lval>=lerr)
			{
				fval = int(val/pow(10.0,lerr)+sign*0.5)*pow(10.0,lerr);
				//count number of trailing zeros
				int digit = 0, digitcount=0, div=1;
				while ((digit==0) && (digitcount<lval))
				{
					digitcount++;
                    div *=10;
					digit = int(fval/pow(10.0,lerr)) % div;
				}
				digitcount --;

				ferr = int(err/pow(10.0,lerr)+0.5);
				if ((lerr<0) && (ferr>9))
				{
					lerr ++;
    				fval = int(fval/pow(10.0,lerr)+sign*0.5)*pow(10.0,lerr);
     				ferr = int(err/pow(10.0,lerr)+0.5);
				}
                
				cout << fval;
				for (i=0; i<digitcount; i++) cout << "0";
				cout << "(" << ferr << ")\n";

				fout << fval;
				for (i=0; i<digitcount; i++) fout << "0";
				fout << "(" << ferr << ")";
			}
			else
			{
				fval = int(val/pow(10.0,lval)+sign*0.5)*pow(10.0,lval);
				ferr = int(err/pow(10.0,lval)+0.5)*pow(10.0,lval);
                cout << fval <<"(" <<ferr<< ")\n";
                fout << fval <<"(" <<ferr<< ")";
			}

		}
	}
    cout << "\n formatted output written to file " << fn_out << '\n';
    cout.flush();
	fin.close();
    fout << " \\\\ \n \\hline \\hline \n";
	fout << "\\end{tabular}\n";
	fout.close();
}
