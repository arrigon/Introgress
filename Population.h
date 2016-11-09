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
    int stock_gamete = -9;

    // Populations
    std::vector <Indiv> indivsPop1;
    std::vector <Indiv> indivsPop2;
    std::vector <int> nUsedGametes;

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
    void loadParams(int nitr,
                    int nphns,
                    double eps,
                    int n_gams,
                    int nind);


    // Initiate parental / hybrid population. TODO: use pointers for that.
    void populateParental(Indiv indinit);                   // with one founder individual

    void populateHybrid(std::vector <Indiv> indsPop1,             // with two parental populations
                        std::vector <Indiv> indsPop2);            // those pops are declared
                                                            // using vectors of individuals


    // Individual related functions
    void getFitnesses(int verbose = 0);       // compute fitnesses of all individuals
    Indiv getOffspring(int verbose = 0);      // produce offsprings, us
    std::vector <Indiv> getAllIndivs();       // output all indivs of pop, as vector


    // Simulation related functions
    double runGenerations(int n_generations,             // Run population for *n_generations*, considering:
                          double mu,                     // mutation rate
                          arma::rowvec phenopt,          // optimal phenotype
                          double omg,                    // selection pressure (omega)
                          double self,                   // selfing rate
                          double bckrte,                 // backcrossing rate (hybrid pops only)
                          std::string outputPath,        // path to outfolder
                          std::string outputPrefix_file, // file prefix where to save genotypes and gamete counts
                          int saveGenotypes = 0,         // save genotypes
                          int save_net_log = 1,          // save genotypes each n generations
                          int verbose = 0);              // print progress to console

    arma::rowvec leftGametes(int verbose = 0);    // count number of gametes remaining after all offsrpings were produced

    // I/O related functions
    void loadNetworks(std::string networks_file);   // Load population from outfile (saved as networks of indivs)
    void savePhenotypes(std::string phens_file);    // Save phenopyte to outfile
    void saveNetworks(std::string networks_file);   // Save networks to outfile
    void saveGameteCounts(std::string gametes_file);     // Save gamete counts to outfile
};


// This is the end of the header guard
#endif
