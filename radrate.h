// $Id: radrate.h 168 2013-04-23 10:13:25Z iamp $
// Stefan Schippers                                 November 2000     
/**
 * \class RADRATE
 * \brief Calculation of hydrogenic dipole transition probabilities 
 * and derived quantities 
 *
 * such as 
 * lifetimes and decay probabilities for all n,l-states up to n=nmax.\n     
 * The n,l specific decay probabilities are calculated from the ion     
 * velocity, vel, and the distances x1 and x2 from the entrance and from
 * the exit of the electron cooler, respectively, to the entrance of  
 * the field ionizing magnet.                                           
 */

class RADRATE
  {
    /**
     * dynamical arrays with size depending on nmax for storage of calculated quantities
     */
    double *transrate_array, *branch_array, *lifetime_array, *pdecay_array;
  int maxn;
  
  public:
  /**
   * constructor calculating hydrogenic transtion rates, lifetimes 
   * and decay probabilities\n
   * for the latter see equation A.1 of Schippers et al., 
   * Astrophys. J. 555 (2001) 1027\n
   * the calculated quantities are stored internally in dynamical arrays
   * @param nmax maximum n to be used in the calculation
   * @param    z effective nuclear charge
   * @param  vel ion velocity  (cm/s)
   * @param   x1 distance from beginnin of cooler to field ionization zone (cm)
   * @param   x2 distance from end of cooler to field ionization zone (cm)\n
   */
     RADRATE(int nmax, double z, double vel, double x1, double x2);

    /**
     * destructor, clears internal storage
     */

    ~RADRATE();

     /** 
      * returns the hydrogenic rate for a n1,l1 -> n2,l2 dipole transtion 
      */
    double trans(int n1, int l1, int n2, int l2) const;

     /** 
      * returns the hydrogenic branching ratio for a n1,l1 -> n2,l2 dipole transtion 
      */
    double branch(int n1, int l1, int n2, int l2) const;

    /** 
     * returns the hydrogenic lifetime of an n,l level
     */
    double life(int n1, int l1) const;

    /**
     * returns the decay probability of an n,l level\n
     * see equation A.1 of Schippers et al., Astrophys. J. 555 (2001) 1027
     */
    double pdecay(int n1, int l1) const ;

    /**
     * destructor, clears internal storage
     */
    void print() const;

  };

