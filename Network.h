// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef NETWORK_H
#define NETWORK_H

#include <iostream>
#include <armadillo>

// This is the content of the .h file, which is where the declarations go
class Network
  /*
    ## Network object ##
     Network is haploid (Wi), endures mutations at a given probability and checks for convergence
     / data containers:
     --- Params
     - w            = gene-to-gene interaction matrix
     - niter_conv   = Nr. of developmental iterations (total)
     - n_phens      = Nr. of iterations to consider for assessing convergence
     - epsilon      = Minimum equilibrium value needed to declare the focal network as convergent
     - mut_rate     = Mutation rate

     --- Variables
     - phen_last    = Stable phenotype (to be accessed via getPhen)
     - conv_test    = Convergence test results (0 = non-convergent)

     / functions:
     - loadParams  = Load user-defined parameters, w is passed via a pointer
     - print       = Shows infos about network
     - getConvg    = Shows convergence test result
     - getPhen     = Shows stable phenotype (if exists, NAN otherwise)
     - getNetwork  = Shows stable phenotype (if exists, NAN otherwise)
     - mutate      = Applies mutations and update network in object. WARNING: we loose original network by doing so.
     - convergence = Applies convergence check
    */

{
private:
    arma::mat w;
    int niter_conv;
    int n_phens;
    int n_genes;
    double epsilon;
    double mut_rate;
    double conv_val;
    double networkFillness;
    arma::rowvec phen_last;
    int conv_test = -9;

public:
 // Constructor definition
    //Network(arma::mat w_ptr, int nitr, int nphns, double eps, double mu);
    void loadParams(arma::mat w_ptr, int nitr, int nphns, double eps, double mu);
    void print();
    int getConvg();
    int rerunConvg();
    arma::rowvec getPhen();
    arma::mat getNetwork();
    void mutate();
    void convergence();
    arma::mat getConvergentNetwork(double netFill);
};
// This is the end of the header guard
#endif
