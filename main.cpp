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


int main()
{
    // initiate variables

    // 16 genes: mu = 0.1, omega = 0.8 / 0.5, full matrix, 100iter conv, eps = 1e-3
    arma::arma_rng::set_seed_random(); // TO UNCOMMENT BEFORE PROD change the random seed
    int n_genes = 49;                    // number of genes in network
    int niter_conv = 100;               // number of iterations in convergence test
    int n_phens = 25;
    int n_gams = 500;
    int n_indiv = 100;
    int n_generations1 = 10;
    int n_generations2 = 10;
    int n_generations3 = 10;

    double epsilon = 1e-3;
    double self_rate = 0;
    double backcross_rate = 0;
    double mut_rate = 1;
    double omega = 0.8;

    rowvec phen_opt = randu<rowvec>(n_genes);


    //// Initiating pop
    // Finding initial convergent network
    mat w(n_genes, n_genes);
    Network net;
    net.loadParams(w,
                    niter_conv,
                    n_phens,
                    epsilon,
                    mut_rate);
    w = net.getConvergentNetwork();



    // Load it into founder individual
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
    cout << "main: Indiv phen_last = " << ind_init.getPhen() << endl;
    Indiv ind_scnd = ind_init;


    // Initiate parental pop with that guy
    Population pop1;
    pop1.loadParams(niter_conv,
                    n_phens,
                    epsilon,
                    n_gams,
                    n_indiv);
    pop1.populateParental(ind_init);



    // Canalize on first phenotype
    phen_opt.load("phen_opt_alien.txt");         // Load phen_opt form outfile

    pop1.runGenerations(n_generations1,          // Run for n_generations1
                        mut_rate,                // - applying mut_rate
                        phen_opt,                // - select on phen_opt
                        omega,                   // - with omega as selection intensity
                        self_rate,               // - selfing rate
                        backcross_rate,          // - backcrossing rate
                        "pop1_dist_to_opt.txt",  // - report distances to phen_opt in outfile
                        1);                      // - verbose mode

    pop1.savePhenotypes("pop0.phens.txt");       // save phenotypes




    // Split in two sub pops
    Population pop2 = pop1;



    // Further evolve Pop1 on phen 1
    phen_opt.load("phen_opt_todd.txt");          // Load phen_opt form outfile

    pop1.runGenerations(n_generations2,          // Run for n_generations2
                        mut_rate,                // - applying mut_rate
                        phen_opt,                // - select on phen_opt
                        omega,                   // - with omega as selection intensity
                        self_rate,               // - selfing rate
                        backcross_rate,          // - backcrossing rate
                        "pop1_dist_to_opt.txt",  // - report distances to phen_opt in outfile
                        1);                      // - verbose mode

    pop1.savePhenotypes("pop1.phens.txt");       // save phenotypes




    // Evolve Pop2 on phen 2
    phen_opt.load("phen_opt_ghost.txt");          // Load phen_opt form outfile

    pop2.runGenerations(n_generations2,          // Run for n_generations2
                        mut_rate,                // - applying mut_rate
                        phen_opt,                // - select on phen_opt
                        omega,                   // - with omega as selection intensity
                        self_rate,               // - selfing rate
                        backcross_rate,          // - backcrossing rate
                        "pop2_dist_to_opt.txt",  // - report distances to phen_opt in outfile
                        1);                      // - verbose mode

    pop2.savePhenotypes("pop2.phens.txt");       // save phenotypes




    //// HYBRIDIZATION PHASE
    // load hybrid population
    Population pop3 = pop1;
    pop3.populateHybrid(pop1.getAllIndivs(),
                        pop2.getAllIndivs());

    // Evolve on Phen 2
    phen_opt.load("phen_opt_ghost.txt");          // Load phen_opt form outfile

    pop3.runGenerations(n_generations3,          // Run for n_generations3
                        mut_rate = 0,                // - applying mut_rate
                        phen_opt,                // - select on phen_opt
                        omega,                   // - with omega as selection intensity
                        self_rate,               // - selfing rate
                        backcross_rate = .3,          // - backcrossing rate
                        "pop3_dist_to_opt.txt",  // - report distances to phen_opt in outfile
                        1);                      // - verbose mode

    pop3.savePhenotypes("pop3.phens.txt");       // save phenotypes




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







//
//    // populate network matrix
//    mat w = randn<mat>(n_genes,n_genes); // initiate network matrix
//
//    // store to Network class
//    Network net1;                      // initiate Network class
//    net1.loadParams(w,
//                    niter_conv,
//                    n_phens,
//                    epsilon,
//                    mut_rate);
//
//    // play around
//    net1.print();
//    net1.convergence();
//    cout << "main: original network = \n" << net1.getNetwork() << endl;   // Get original network                                               // apply convergence test
//    cout << "main: conv_test 1 = " << net1.getConvg() << endl;            // Get convergence test 1
//    cout << "main: phen_last 1 = " << net1.getPhen() << endl;             // Get phenotype 1
//
//
//    // mutate and get new network out of it
//    net1.mutate();
//    mat w2 = net1.getNetwork();
//
//    Network net2;                      // initiate Network class
//    net2.loadParams(w2,
//                    niter_conv,
//                    n_phens,
//                    epsilon,
//                    mut_rate);
//    net2.print();
//    net2.convergence();
//    cout << "main: mutated network = \n" << net2.getNetwork() << endl;    // Get network
//    cout << "main: conv_test 2 = " << net2.getConvg() << endl;            // Get convergence test 2
//    cout << "main: phen_last 2 = " << net2.getPhen() << endl;             // Get phenotype 2
//
//
//    // Playing with pointers
//    Network *net2_ptr = &net2;  // store net2 as a pointer
//    net2_ptr->print();          // access to Network functions
//    cout << "main: network 2, via ptr = \n" << net2_ptr->getNetwork() << endl;    // Get network
//    net2_ptr->mutate();
//    mat net3 = net2_ptr->getNetwork();                                     // Mutate network
//    cout << "main: network 3, mutated via net2_ptr = \n" << net3 << endl;         // Get network
//
//
//    // Playing with dynamic allocation
//    int nstable = 0;
//    int attempts = 0;
//
//    while(nstable < 10)
//    {
//        net1.mutate();
//        mat w_mut = net1.getNetwork();
//        Network *net_ptr = new Network();                      // initiate Network class, with new pointer
//        net_ptr->loadParams(w_mut,
//                            niter_conv,
//                            n_phens,
//                            epsilon,
//                            mut_rate);
//        int conv_test = net_ptr->convergence();
//        if(conv_test == 1)
//        {
//            cout << "attmpt: " << attempts << " " << nstable << " " << net_ptr << endl;
//            nstable++;
//
//        }
//        else
//        {
//            delete net_ptr;
//        }
//        attempts++;
//    }
//
//
//    // Playing mith Indiv object
//    Indiv ind1;                      // initiate Indiv class
//    ind1.rdmStart(n_genes,
//                  niter_conv,        // w1 and w2 are declared internally, to use for debug purposes only
//                  n_phens,
//                  epsilon,
//                  mut_rate,
//                  n_gams,
//                  phen_opt,
//                  omega);
//
//    // print gene networks
//    cout << "main: w1 = \n" << ind1.getNetwork1() << endl;
//    cout << "main: w2 = \n" << ind1.getNetwork2() << endl;
//
//    // get Indiv phenotype and fitness
//    rowvec w_phen = ind1.getPhen();             // Get phenotype
//    cout << "main: phen_opt = \n" << phen_opt << endl;
//    cout << "main: phen = \n" << w_phen << endl;
//
//    double ind1_fitn = ind1.getFitn();          // Get fitness
//    cout << "main: ind1_fitn = " << ind1_fitn << endl;
//
//    // produce gamete
//    mat w_gamete = ind1.getGamete();                  // get gamete
//    cout << "main: w_gamete = \n" << w_gamete << endl;
//    cout << "main: ind1 remaining gametes = " << ind1.getGamCnt() << endl;


//// save to outfile
//    string test = "pop1.phens.txt";
//    pop1.savePhenotypes(test);
//    pop1.saveNetworks("pop1.nets.txt");
//
//
//    // Load population from existing file
//    Population pop2;
//    Indiv ind_null;
//    pop2.loadParams(ind_null,
//                    ind_null,
//                    niter_conv,
//                    n_phens,
//                    epsilon,
//                    mut_rate,
//                    n_gams,
//                    phen_opt,
//                    omega,
//                    n_indiv,
//                    self_rate,
//                    backcross_rate);
//    pop2.loadNetworks("pop1.nets.txt");
//    pop2.saveNetworks("pop2.nets.txt");
//    pop2.getFitnesses(0);
//    pop2.savePhenotypes("pop2.phens.txt");
//    pop2.runGenerations(1);
//
