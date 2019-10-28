// $Id: readxsec.h 497 2018-02-02 17:40:30Z iamp $
void info_DRtheory(void);

int open_DRtheory(char *filename, int &theo_mode, int &number_of_levels);

void read_DRtheory(char *filename, int theo_mode, int level_number, double *energy, double *strength, 
                 double *wLorentz, int npts, double &emin, double &emax);


int open_peak_file(char *filename, int &theo_mode, int &number_of_levels);

int read_peak_file(char *filename, int theo_mode, int level_number, double *energy, double *strength, 
		    double *qFano, double *wLorentz, double *wGauss, int npts, double &emin, double &emax);
