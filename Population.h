// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef POPULATION_H
#define POPULATION_H

#include <armadillo>
#include <random>
#include <vector>


class Population
{
    /*
    ## Population object ##
    / Population:
      - is made of individuals, via pointers to Indiv objects
      - ranks those individuals according to their fitness
      - produces offsprings, by pairing Indivs according to fitness and/or breeding instructions

    / data containers
      - stock_indivs = 2 columns array
                       - pointers to Indivs
                       - associated fitness
                       - pop of origin (for backcrossing)

      - scnd_parent  = forced cross, via a pointer to the second parent

    / functions:
      - loadParams    = Loads user-defined parameters
      - populate      = Fills population with ind_init (one species mode)
      - backcross     = Plays with two species
      - getOffspring  = function producing offsprings on demand,
                        used to populate the next generation
                        this function updates stock_indivs if an individual exhausts its gamete stock
                        (except for the scnd_parent).
    */
private:
    // Network related
    arma::mat w1, w2, w_gamete;

    // Individuals
    std::vector <Indiv> myIndivs;
    std::vector <int> myNames;
    std::vector <arma::mat> myGametes;
    std::vector <double> myFitnesses;
    Indiv ind_init;
    Indiv ind_scndprt;
    int stock_gamete = -9;

    // General params
    // -> pop params
    int n_indiv;
    int n_generations;
    double backcross_rate; // backcross rate

    // -> convergence params
    int niter_conv;
    int n_phens;
    double epsilon;
    arma::rowvec phen_last;
    int conv_test = -9;
    double mut_rate;
    double self_rate;

    // -> fitness related
    arma::rowvec phen_opt;
    double omega;
    double fitvalue = 0;

    // Functions
public:
    // Constructor definition
    void loadParams(Indiv indinit,
                    Indiv indscnd,
                    int nitr,
                    int nphns,
                    double eps,
                    double mu,
                    int n_gams,
                    arma::rowvec phenopt,
                    double omg,
                    int nind,
                    double self,
                    double bckrte);

    void newPhenOpt(arma::rowvec phenopt);
    void populateIni();
    void spanGametes();
    void getFitnesses(int verbose = 0);
    Indiv getOffspring(int verbose = 0);
    void runGenerations(int n_generations, int verbose = 0);
    void savePhenotypes(std::string phens_file);
    void saveNetworks(std::string networks_file);
    void loadNetworks(std::string networks_file);
};


// This is the end of the header guard
#endif
