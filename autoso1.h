// $Id: autoso1.h 168 2013-04-23 10:13:25Z iamp $
// Author: Stefan.E.Schippers@iamp.physik.uni-giessen.de 2006-07-18

/** \class AUTOSO1
* \brief Reads, holds and manipulates a nl-block 
*  of data from an autostructure *.o1 file
*
*/
#include <fstream>

using namespace std;

void TestAutosO1(void);

class AUTOSO1
{
    int nAA, nLV, nAR;

    int *AA_CF1, *AA_LV1, *AA_W1, *AA_CF2, *AA_LV2;
    double *AA_rate, *AA_energy;

    int *LV_K, *LV_LV, *LV_T, *LV_S, *LV_L, *LV_J, *LV_CF;
    double *LV_energy;  

    int *AR_CF1, *AR_LV1, *AR_W1, *AR_CF2, *AR_LV2, *AR_W2;
    double *AR_rate, *AR_energy;

    int nRydberg, lRydberg;
    double Zeff;

    fstream fo1;
  public:

     ~AUTOSO1();

     int open_file(char* filename); /**< opens the *.o1 file */

     void close_file(void);         /**< closes the *.o1 file */

     int read_nl(int nRyd, int lRyd); /**< reads a block of data 
                            belonging to a given nl-Rydberg state */ 

     int summed_rates(int nminHydro, int CF1); /**< manipulates data
   			    beloging to a given nl-Rydberg state */ 
  };






