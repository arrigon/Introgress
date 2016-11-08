#include "Maths.h"
#include "Network.h"
#include "Indiv.h"
#include "Population.h"
//#include <iostream>
//#include <armadillo>
//#include <random>
//#include <vector>

using namespace std;
using namespace arma;


int main(int argc, char *argv[])
{
    arma::arma_rng::set_seed_random(); // TO UNCOMMENT BEFORE PROD change the random seed

    //DEBUG SWITCH ///////
    int debugMode = 1;  // WARNING: debugMode = 0 for PROD !!!
    //////////////////////


    //// Initiate variables
    int n_genes, niter_conv, n_phens, n_gams, n_indiv, n_generations1, n_generations2, n_generations3, connectGradient;                 // number of genes in network
    double networkFillness, epsilon, self_rate, backcross_rate, mut_rate, omega, mindist;      // proportion of non-null interactions in the network


    // initiate variables DEBUG MODE ONLY
    if(debugMode == 1)
    {
        n_genes = 49;                 // number of genes in network
        networkFillness = .5;      // proportion of non-null interactions in the network
        connectGradient = 1;          // set a connectivity gradient in the gene network

        niter_conv = 100;             // nb of iterations in convergence test
        n_phens = 25;                 // nb of phenotypes we keep in convergence tests
        epsilon = 1e-5;            // convergence criterion

        n_gams = 500;                 // nb of gametes available per individual
        n_indiv = 100;                // nb of individuals per population

        n_generations1 = 50;          // nb of generations in training phase
        n_generations2 = 200;         // nb of generations in divergence phase
        n_generations3 = 30;          // nb of generations in hybridization phase

        self_rate = 0;             // selfing rate
        backcross_rate = 0;        // backcrossing rate
        mut_rate = 0.05;           // mutation rate (during phases 1 and 2, set to 0 via hard coding in hybridization phase)
        omega = 0.8;               // selection coefficient
        mindist = 0.05;            // (used instead of n_generations1, = minimum phenotypic distance to reach before moving to divergence phase)
    }


    // Initiate variables PRODUCTION MODE
    if(debugMode == 0)
    {
        n_genes = atoi(argv[1]);                 // number of genes in network
        networkFillness = atof(argv[2]);      // proportion of non-null interactions in the network
        connectGradient = 1;                     // set a connectivity gradient in the gene network

        niter_conv = atoi(argv[3]);             // nb of iterations in convergence test
        n_phens = atoi(argv[4]);                // nb of phenotypes we keep in convergence tests
        epsilon = atof(argv[10]);            // convergence criterion

        n_gams = atoi(argv[5]);                 // nb of gametes available per individual
        n_indiv = atoi(argv[6]);                // nb of individuals per population

        n_generations1 = atoi(argv[7]);          // nb of generations in training phase
        n_generations2 = 200;         // nb of generations in divergence phase
        n_generations3 = 30;          // nb of generations in hybridization phase

        self_rate = atof(argv[11]);             // selfing rate
        backcross_rate = 0;        // backcrossing rate
        mut_rate = atof(argv[13]);  // 13 ou 14 ??           // mutation rate (during phases 1 and 2, set to 0 via hard coding in hybridization phase)
        omega = atof(argv[14]);               // selection coefficient
        mindist = 0.05;            // (used instead of n_generations1, = minimum phenotypic distance to reach before moving to divergence phase)
    }


    //// Initiating pop
    // Finding initial convergent network
    cout << "main: n_genes =\n" << n_genes << endl;
    mat w(n_genes, n_genes);
    Network net;
    net.loadParams(w,
                    niter_conv,
                    n_phens,
                    epsilon,
                    mut_rate);
    w = net.getConvergentNetwork(networkFillness, connectGradient);
    cout << "main: convergent network =\n" << w << endl;



//    // Load it into founder individual
//    Indiv ind_init;                      // initiate Indiv class
//    ind_init.loadParams(w,
//                        w,
//                        niter_conv,        // w1 and w2 are declared internally, to use for debug purposes only
//                        n_phens,
//                        epsilon,
//                        mut_rate,
//                        n_gams,
//                        phen_opt,
//                        omega);
//    cout << "main: Indiv phen_last = " << ind_init.getPhen() << endl;
//    Indiv ind_scnd = ind_init;
//
//
//    // Initiate parental pop with that guy
//    Population pop1;
//    pop1.loadParams(niter_conv,
//                    n_phens,
//                    epsilon,
//                    n_gams,
//                    n_indiv);
//    pop1.populateParental(ind_init);
//
//
//
//    // Canalize on first phenotype
//    rowvec phen_opt(n_genes);
//    phen_opt.load("phen_opt_alien.txt");         // Load phen_opt form outfile
//
//    double distToOpt = 1;
//    while (distToOpt > mindist)
//    {
//        distToOpt = pop1.runGenerations(1,                          // Run for 1 generation
//                                        mut_rate,                   // - applying mut_rate
//                                        phen_opt,                   // - select on phen_opt
//                                        omega,                      // - with omega as selection intensity
//                                        self_rate,                  // - selfing rate
//                                        backcross_rate,             // - backcrossing rate
//                                        "pop0",           // - genotypes outfile
//                                        1,                          // - save genotypes to outfile (0 = F, 1 = T)
//                                        1);                         // - verbose mode
//    }


//    pop1.savePhenotypes("pop0.phens.txt");       // save phenotypes
//
//
//
//
//    // Split in two sub pops
//    Population pop2 = pop1;
//
//
//
//    // Further evolve Pop1 on phen 1
//    phen_opt.load("phen_opt_todd.txt");          // Load phen_opt form outfile
//
//    pop1.runGenerations(n_generations2,          // Run for n_generations2
//                        mut_rate,                // - applying mut_rate
//                        phen_opt,                // - select on phen_opt
//                        omega,                   // - with omega as selection intensity
//                        self_rate,               // - selfing rate
//                        backcross_rate,          // - backcrossing rate
//                        "pop1_genotypes.txt",    // - report distances to phen_opt in outfile
//                        0,                       // - save genotypes to outfile (0 = F, 1 = T)
//                        1);                      // - verbose mode
//
//    pop1.savePhenotypes("pop1.phens.txt");       // save phenotypes
//
//
//
//
//    // Evolve Pop2 on phen 2
//    phen_opt.load("phen_opt_ghost.txt");          // Load phen_opt form outfile
//
//    pop2.runGenerations(n_generations2,          // Run for n_generations2
//                        mut_rate,                // - applying mut_rate
//                        phen_opt,                // - select on phen_opt
//                        omega,                   // - with omega as selection intensity
//                        self_rate,               // - selfing rate
//                        backcross_rate,          // - backcrossing rate
//                        "pop2_genotypes.txt",  // - report distances to phen_opt in outfile
//                        0,                       // - save genotypes to outfile (0 = F, 1 = T)
//                        1);                      // - verbose mode
//
//    pop2.savePhenotypes("pop2.phens.txt");       // save phenotypes
//
//
//
//    //// HYBRIDIZATION PHASE
//    // load hybrid population
//    Population pop3;
//    pop3.loadParams(niter_conv,
//                    n_phens,
//                    epsilon,
//                    n_gams,
//                    n_indiv);
//    pop3.populateHybrid(pop1.getAllIndivs(),
//                        pop2.getAllIndivs());
//    pop3.saveNetworks("pop3_F1genotypes.txt");
//
//    // Evolve on Phen 1, and backcross on Pop2
//    phen_opt.load("phen_opt_todd.txt");          // Load phen_opt form outfile
//
//    pop3.runGenerations(n_generations3,          // Run for n_generations3
//                        mut_rate = 0,                // - applying mut_rate
//                        phen_opt,                // - select on phen_opt
//                        omega,                   // - with omega as selection intensity
//                        self_rate,               // - selfing rate
//                        backcross_rate = 0.95,     // - backcrossing rate, we backcross on **Pop2**
//                        "pop3_dist_to_opt.txt",  // - report distances to phen_opt in outfile
//                        1);                      // - verbose mode
//
//    pop3.savePhenotypes("pop3.phens.txt");       // save phenotypes
//    pop3.saveNetworks("pop3_BCgenotypes.txt");
//
//

//    cout << "main: Indiv W_last = \n" << ind_final.getNetwork1() << endl;
//
//    rowvec phen_trained = ind_final.getPhen();
//    phen_trained.save("phen_trained1.txt");
//    mat w_trained = ind_final.getNetwork1();
//    w_trained.save("Net_trained1.txt");

     //cout << "Final fitnesses: " << endl;
     //pop1.getFitnesses(1);

    // end of main
    cout << "main: end OK OK OK OK OK \n";

    return 0;
}

