// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef INDIV_H
#define INDIV_H

#include <armadillo>
#include <random>

class Indiv
{
    /*
    ## Indiv object ##
     Indiv:
     - is diploid genotype (via pointers to two network objects, w1 and w2),
     - shows an associated phenotype (getPhen) and fitness (getFitn),
     - applies meiosis through the following steps:
       --> 1. recombine parental networks
       --> 2. apply mutations (according to mut_rate, via network.mutate method)
       --> 3. checks for gamete viablity (via network.getConvg method)
       --> 4. return gamete if viable, w.fill(NAN) if not
       --> 5. decrement stock of available gametes

     / data containers:
     - w1, w2        = pair of networks objects
     - stock_gamete  = number of available gametes
     - phen_last     = phenotype
     - fitn_last     = fitness value

     / functions
     - rdmStart       = Initiates random individual (debug only)
     - loadParams     = Loads user-defined parameters
     - getNetwork1|2  = Get specimen network w1 or w2
     - getPhen        = Get specimen phenotype
     - getFitn        = Get specimen fitness
     - getGamete      = Produce one gamete out of available pool, applies meiosis
                        -> takes "infinite = 0" or infinite = 1 param, to specifiy if the focal individual
                           is an infinite provider of gametes
                        -> return either gamete network or matrix filled with NAN if non convergent
    */


private:
    arma::mat w1, w2, w_gamete;
    // Network net_ind_tmp;
    Network net_gamete;
    int niter_conv;
    int n_phens;
    double epsilon;
    double mut_rate;
    arma::rowvec phen_last;
    arma::rowvec phen_opt;
    double fitvalue;
    int conv_test = -9;
    int stock_gamete;
    double omega;

    // Functions
public:
    void rdmStart(int n_genes,
                  int nitr,
                  int nphns,
                  double eps,
                  double mu,
                  int n_gams,
                  arma::rowvec phenopt,
                  double omg);

    void loadParams(arma::mat w_1,
                    arma::mat w_2,
                    int nitr,
                    int nphns,
                    double eps,
                    double mu,
                    int n_gams,
                    arma::rowvec phenopt,
                    double omg);
    arma::mat getNetwork1();
    arma::mat getNetwork2();
    int getGameteCnt();
    arma::rowvec getPhen();
    double getFitn();
    arma::mat getGamete(int infinite);
};

// This is the end of the header guard
#endif
