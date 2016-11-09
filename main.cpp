#include "Maths.h"
#include "Network.h"
#include "Indiv.h"
#include "Population.h"
// #include <boost/filesystem.hpp>


//#include <iostream>
//#include <armadillo>
//#include <random>
//#include <vector>

using namespace std;
using namespace arma;


int main(int argc, char *argv[])
{
    arma::arma_rng::set_seed_random(); // TO UNCOMMENT BEFORE PROD change the random seed

    //DEBUG SWITCHES ///////
    int debugMode = 1;     // WARNING: debugMode = 0 for PROD !!!
    int verboseMain = 0;             // if main must be verbose or not
    int verboseRunGenerations = 0;   // if main must be verbose or not
    ////////////////////////


    //// Initiate variables
    int n_genes, niter_conv, n_phens, n_gams, n_indiv, n_generations1, n_generations2, n_generations3, connectGradient, netlog;                 // number of genes in network
    double networkFillness, epsilon, self_rate, backcross_rate, mut_rate, omega, mindist;      // proportion of non-null interactions in the network
    string popname, outputPath;

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

        n_generations1 = 10;          // nb of generations in training phase
        n_generations2 = 10;         // nb of generations in divergence phase
        n_generations3 = 10;          // nb of generations in hybridization phase

        self_rate = 0;                // selfing rate
        backcross_rate = 0;           // backcrossing rate
        mut_rate = 0.05;              // mutation rate (during phases 1 and 2, set to 0 via hard coding in hybridization phase)
        omega = 0.8;                  // selection coefficient
        mindist = 0.05;               // (used instead of n_generations1, = minimum phenotypic distance to reach before moving to divergence phase)
        outputPath = "testOutputs/";  // Folder where to save outputs
        netlog = 10;                  // Interval of generations when to save outputs
    }


    // Initiate variables PRODUCTION MODE
    if(debugMode == 0)
    {
        n_genes = atoi(argv[1]);                // number of genes in network
        networkFillness = atof(argv[2]);        // proportion of non-null interactions in the network
        connectGradient = 1;                    // set a connectivity gradient in the gene network

        niter_conv = atoi(argv[3]);             // nb of iterations in convergence test
        n_phens = atoi(argv[4]);                // nb of phenotypes we keep in convergence tests
        epsilon = atof(argv[10]);               // convergence criterion

        n_gams = atoi(argv[5]);                 // nb of gametes available per individual
        n_indiv = atoi(argv[6]);                // nb of individuals per population

        n_generations1 = atoi(argv[7]);         // nb of generations in training phase
        n_generations2 = atoi(argv[8]);         // nb of generations in divergence phase
        n_generations3 = atoi(argv[9]);         // nb of generations in hybridization phase

        self_rate = atof(argv[11]);             // selfing rate
        backcross_rate = atof(argv[12]);        // backcrossing rate
        mut_rate = atof(argv[13]);              // mutation rate (during phases 1 and 2, set to 0 via hard coding in hybridization phase)
        omega = atof(argv[14]);                 // selection coefficient
        mindist = 0.05;                         // minimum phenotypic distance to reach before moving to divergence phase (used instead of n_generations1)
        outputPath = argv[15];                        // Folder where to save outputs
        netlog = atoi(argv[16]);                // Interval of generations when to save outputs
    }


    //// Preparing outfolder
//    boost::filesystem::path dir(outputPath);
//
//    if(!(boost::filesystem::exists(dir))){
//        if(verboseMain == 1) cout << "Doesn't Exists" << endl;
//
//        if (boost::filesystem::create_directory(dir))
//        {
//        if(verboseMain == 1) cout << "main: creating " << outputPath << endl;
//        }
//    }





    ///////////////////////
    //// CANALIZATION PHASE
    ///////////////////////
    cout << "main: *** CANALIZATION PHASE ***" << endl;

    // Finding initial convergent network
    mat w(n_genes, n_genes);
    Network net;
    net.loadParams(w,
                    niter_conv,
                    n_phens,
                    epsilon,
                    mut_rate);
    w = net.getConvergentNetwork(networkFillness, connectGradient);
    if(verboseMain == 1) cout << "main: initial convergent network =\n" << w << endl;


    // Load it into founder individual
    rowvec phen_opt = randu<rowvec>(n_genes);
    Indiv ind_init;                      // initiate Indiv class
    ind_init.loadParams(w,
                        w,
                        niter_conv,        // w1 and w2 are declared internally, to use for debug purposes only
                        n_phens,
                        epsilon,
                        mut_rate,
                        n_gams,
                        phen_opt,
                        omega);

    // Get phenotype of that guy, declare it as optimal phenotype
    rowvec phen_opt0 = ind_init.getPhen();
    if(verboseMain == 1) cout << "main: Indiv phen_last = " << ind_init.getPhen() << endl;


    // Initiate parental pop with that guy
    Population pop1;
    pop1.loadParams(niter_conv,
                    n_phens,
                    epsilon,
                    n_gams,
                    n_indiv);
    pop1.populateParental(ind_init);


    // Canalize on this initial phenotype, for n_generations1
    // phen_opt0.load("phen_opt_alien.txt");        // IF WILLING TO USE PRECOMPUTED PHENOTYPES
    pop1.runGenerations(n_generations1,             // Run for n_generations1
                        mut_rate,                   // - applying mut_rate
                        phen_opt0,                  // - select on phen_opt0,
                        omega,                      // - with omega as selection intensity
                        self_rate,                  // - selfing rate
                        0,                          // - NO BACKCROSSING IN PARENTAL POPS
                        outputPath,                 // - outfiles path
                        "pop0",                     // - outfiles prefix
                        0,                          // - save genotypes to outfile (0 = F, 1 = T)
                        netlog,                     // - save every netlog generations
                        verboseRunGenerations);     // - verbose mode





    ///////////////////////
    //// DIVERGENCE PHASE
    ///////////////////////
    cout << "main: *** DIVERGENCE PHASE ***" << endl;

    // Split in two sub pops
    Population pop2 = pop1;


    // POP1
    cout << "main: ---> Pop1 " << endl;

    // Evolve Pop1 on phen 1 (random)
    rowvec phen_opt1 = randu<rowvec>(n_genes);
    // phen_opt1.load("phen_opt_ghost.txt");        // IF WILLING TO USE PRECOMPUTED PHENOTYPES

    // Evolve Pop1 on phen_opt1, for n_generations2
    pop1.runGenerations(n_generations2,             // Run for n_generations2
                        mut_rate,                   // - applying mut_rate
                        phen_opt1,                  // - select on phen_opt0,
                        omega,                      // - with omega as selection intensity
                        self_rate,                  // - selfing rate
                        0,                          // - NO BACKCROSSING IN PARENTAL POPS
                        outputPath,                 // - outfiles path
                        "pop1",                     // - outfiles prefix
                        1,                          // - save genotypes to outfile (0 = F, 1 = T)
                        netlog,                     // - save every netlog generations
                        verboseRunGenerations);     // - verbose mode



    // POP2
    cout << "main: ---> Pop2 " << endl;

    // Evolve Pop2 on phen 2 (random)
    rowvec phen_opt2 = randu<rowvec>(n_genes);
    // phen_opt2.load("phen_opt_todd.txt");        // IF WILLING TO USE PRECOMPUTED PHENOTYPES

    // Evolve Pop1 on phen_opt1, for n_generations2
    pop2.runGenerations(n_generations2,             // Run for n_generations2
                        mut_rate,                   // - applying mut_rate
                        phen_opt2,                  // - select on phen_opt0,
                        omega,                      // - with omega as selection intensity
                        self_rate,                  // - selfing rate
                        0,                          // - NO BACKCROSSING IN PARENTAL POPS
                        outputPath,                 // - outfiles path
                        "pop2",                     // - outfiles prefix
                        1,                          // - save genotypes to outfile (0 = F, 1 = T)
                        netlog,                     // - save every netlog generations
                        verboseRunGenerations);     // - verbose mode





    ////////////////////////
    //// HYBRIDIZATION PHASE
    ////////////////////////
    cout << "main: *** HYBRIDIZATION PHASE ***" << endl;

    // load hybrid population
    Population pop3;
    pop3.loadParams(niter_conv,
                    n_phens,
                    epsilon,
                    n_gams,
                    n_indiv);
    pop3.populateHybrid(pop1.getAllIndivs(),
                        pop2.getAllIndivs());

    // Make sure we save the F1 hybrids
    string F1_file = outputPath;
    F1_file = F1_file + "/Genotypes_F1_pop3.txt";
    pop3.saveNetworks(F1_file);


    // Evolve Pop3 on phen_opt1, backcross on Pop2, for n_generations3
    pop3.runGenerations(n_generations3,             // Run for n_generations2
                        mut_rate = 0,               // - NO MUTATIONS IN HYBRID POP
                        phen_opt1,                  // - SELECT ON phen_opt1
                        omega,                      // - with omega as selection intensity
                        self_rate,                  // - selfing rate
                        backcross_rate,             // - backcrossing rate
                        outputPath,                 // - outfiles path
                        "pop3",                     // - outfiles prefix
                        1,                          // - save genotypes to outfile (0 = F, 1 = T)
                        1,                          // - SAVE EVERY GENERATION FOR HYBRID POP
                        verboseRunGenerations);     // - verbose mode





    // end of main
    cout << "main: *** DONE ***";
    return 0;
}

