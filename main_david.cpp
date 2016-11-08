#include "Maths.h"
#include "Network.h"
#include "Indiv.h"
#include "Population.h"
#include <iostream>
#include <armadillo>
#include <random>
#include <vector>

using namespace std;
using namespace arma;

//    int n_genes = 49;                     // number of genes in network;                                                                  argv[1]
//    double networkFillness = .75;         // proportion of non-null interactions in the network;                                          argv[2]
//    int niter_conv = 100;                 // number of iterations in convergence test ;                                                   argv[3]
//    int n_phens = 25;                                                                                                                     argv[4]
//    int n_gams = 500;                                                                                                                     argv[5]
//    int n_indiv = 100;                                                                                                                    argv[6]
//    int n_generations1 = 100;                                                                                                             argv[7]
//    int n_generations2 = 200;                                                                                                             argv[8]
//    int n_generations3 = 100;                                                                                                             argv[9]
//    double epsilon = 1e-3;                                                                                                                argv[10]
//    double self_rate = 0;                                                                                                                 argv[11]
double backcross_rate;                      // still need to be defined because of pre-Backcrossing generation where backcross= 0;          argv[12]
double mut_rate;                            // still need to be defined because of Backcrossing scenatrio where mutation= 0;                argv[13]
//    double omega = 0.8;                                                                                                                   argv[14]
//                                          // dossier pour stocker les fichiers de chaque réplicats                                        argv[15]

int main(int argc,char *argv[])
{
//        std::string string(arg[15]);
                                                    //     initiate variables
                                                    // 16 genes: mu = 0.1, omega = 0.8 / 0.5, full matrix, 100iter conv, eps = 1e-3
        arma::arma_rng::set_seed_random();          // TO UNCOMMENT BEFORE PROD change the random seed

        rowvec phen_opt = randu<rowvec>(atoi(argv[1]));

                                                    //// Initiating pop
                                                    // Finding initial convergent network
        mat w(atoi(argv[1]),atoi(argv[1]));
        Network net;
        net.loadParams(w,
                     atoi(argv[3]),
                     atoi(argv[4]),
                     atof(argv[10]),
                     atof(argv[14]));

        cout << atof(argv[2]) << endl;
        w = net.getConvergentNetwork(atof(argv[2]));


                                                    // Load it into founder individual
        Indiv ind_init;                             // initiate Indiv class
        ind_init.loadParams(w,
                          w,
                          atoi(argv[3]),            // w1 and w2 are declared internally, to use for debug purposes only
                          atoi(argv[4]),
                          atof(argv[10]),
                          atof(argv[13]),
                          atoi(argv[5]),
                          phen_opt,
                          atof(argv[14]));
        cout << "main: Indiv phen_last = " << ind_init.getPhen() << endl;
        Indiv ind_scnd = ind_init;


                                                    // Initiate parental pop with that guy
        Population pop1;
        pop1.loadParams(atoi(argv[3]),
                      atoi(argv[4]),
                      atof(argv[10]),
                      atoi(argv[5]),
                      atoi(argv[6]));
        pop1.populateParental(ind_init);

                                                    // Canalize on first phenotype
        phen_opt.load("phen_opt_alien.txt");        // Load phen_opt form outfile


        std::stringstream save_dist_to_opt0_file;    // allows to save the file in the correct replicate folder
        save_dist_to_opt0_file << argv[15] << "/pop0_dist_to_opt.txt";
        pop1.runGenerations(atoi(argv[7]),          // Run for n_generations1
                          atof(argv[13]),           // - applying mut_rate
                          phen_opt,                 // - select on phen_opt
                          atof(argv[14]),           // - with omega as selection intensity
                          atof(argv[11]),           // - selfing rate
                          backcross_rate = 0,           // - backcrossing rate
                          save_dist_to_opt0_file.str(),   // - report distances to phen_opt in outfile
                          argv[15],                 // File where to save the genotypes
                          "pop",
                          30,
                          1);                       // - verbose mode

        std::stringstream save_phen0_file;          // allows to save phenotypes in the correct replicate folder
        save_phen0_file << argv[15] << "/pop0.phens.txt";
        pop1.savePhenotypes(save_phen0_file.str());      // allows to save phenotypes in the correct replicate folder
        std::stringstream save_gen0_file;
        save_gen0_file << argv[15] << "/pop0.genotypes.txt";
        pop1.saveNetworks(save_gen0_file.str());

        Population pop2 = pop1;                     // Split in two sub pops

                                                    // Further evolve Pop1 on phen 1
        phen_opt.load("phen_opt_todd.txt");         // Load phen_opt form outfile



        std::stringstream save_dist_to_opt1_file;    // allows to save the file in the correct replicate folder
        save_dist_to_opt1_file << argv[15] << "/pop1_dist_to_opt.txt";
        save_phen0_file << argv[15] << "/pop0.phens.txt";
        pop1.runGenerations(atoi(argv[8]),          // Run for n_generations2
                          atof(argv[13]),           // - applying mut_rate
                          phen_opt,                 // - select on phen_opt
                          atof(argv[14]),           // - with omega as selection intensity
                          atof(argv[11]),           // - selfing rate
                          backcross_rate = 0,           // - backcrossing rate
                          save_dist_to_opt1_file.str(),   // - report distances to phen_opt in outfile
                          argv[15],                 // File where to save the genotypes
                          "pop1",
                          30,
                          1);                       // - verbose mode

        std::stringstream save_phen1_file;
        save_phen1_file << argv[15] << "/pop1.phens.txt";
        pop1.savePhenotypes(save_phen1_file.str());      // save phenotypes

                                                    // Evolve Pop2 on phen 2
        phen_opt.load("phen_opt_ghost.txt");        // Load phen_opt form outfile


        std::stringstream save_dist_to_opt2_file;    // allows to save the file in the correct folder
        save_dist_to_opt2_file << argv[15] << "/pop2_dist_to_opt.txt";
        pop2.runGenerations(atoi(argv[8]),          // Run for n_generations2
                          atof(argv[13]),           // - applying mut_rate
                          phen_opt,                 // - select on phen_opt
                          atof(argv[14]),           // - with omega as selection intensity
                          atof(argv[11]),           // - selfing rate
                          backcross_rate = 0,           // - backcrossing rate
                          save_dist_to_opt2_file.str(),   // - report distances to phen_opt in outfile
                          argv[15],                 // File where to save the genotypes
                          "pop2",
                          30,
                          1);                       // - verbose mode
        std::stringstream save_phen2_file;
        save_phen2_file << argv[15] << "/pop2.phens.txt";
        pop2.savePhenotypes(save_phen2_file.str());      // save phenotypes



                                                    //// HYBRIDIZATION PHASE
                                                    // load hybrid population
        Population pop3;
        pop3.loadParams(atoi(argv[3]),
                      atoi(argv[4]),
                      atof(argv[10]),
                      atoi(argv[5]),
                      atoi(argv[6]));
        pop3.populateHybrid(pop1.getAllIndivs(),
                          pop2.getAllIndivs());

        std::stringstream save_genF1_file;
        save_genF1_file << argv[15] << "/pop3_F1genotypes.txt";
        pop3.saveNetworks(save_genF1_file.str());


                                                    // Evolve on Phen 1, and backcross on Pop2
        phen_opt.load("phen_opt_todd.txt");         // Load phen_opt form outfile

        std::stringstream save_dist_to_opt3_file;   // allows to save the file in the correct folder
        save_dist_to_opt3_file << argv[15] << "/pop3_dist_to_opt.txt";
        pop3.runGenerations(atoi(argv[9]),          // Run for n_generations3
                          mut_rate = 0,             // - applying mut_rate
                          phen_opt,                 // - select on phen_opt
                          atof(argv[14]),           // - with omega as selection intensity
                          atof(argv[11]),           // - selfing rate
                          atof(argv[12]),    // - backcrossing rate, we backcross on **Pop2**
                          save_dist_to_opt3_file.str(),   // - report distances to phen_opt in outfile
                          argv[15],
                          "pop3", // File where to save the genotypes
                          1,
                          1);                       // - verbose mode

        std::stringstream save_genBC3_file;
        save_genBC3_file << argv[15] << "/pop3_BCgenotypes.txt";
        pop3.saveNetworks(save_genBC3_file.str());          // save phenotypes

                                                    //    cout << "main: Indiv W_last = \n" << ind_final.getNetwork1() << endl;
                                                    //    rowvec phen_trained = ind_final.getPhen();
                                                    //    phen_trained.save("phen_trained1.txt");
                                                    //    mat w_trained = ind_final.getNetwork1();
                                                    //    w_trained.save("Net_trained1.txt");
                                                    //     cout << "Final fitnesses: " << endl;
        pop1.getFitnesses(1);
        std::stringstream save_phen3_file;
        save_phen3_file << argv[15] << "/pop3.phens.txt";
        pop3.savePhenotypes(save_phen3_file.str()); // save phenotypes

        //    end of main
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
