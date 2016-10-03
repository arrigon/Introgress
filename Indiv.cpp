#include "Network.h"
#include "Indiv.h"
#include "Maths.h"
#include <armadillo>
#include <random>


using namespace std;
using namespace arma;

void Indiv::rdmStart(int n_genes,
                     int nitr,
                     int nphns,
                     double eps,
                     double mu,
                     int n_gams,
                     rowvec phenopt,
                     double omg)
{
    w1 = randn<mat>(n_genes,n_genes);   // gene-to-gene interaction matrix
    w2 = w1;                            // gene-to-gene interaction matrix
    niter_conv = nitr;            // number of iterations to check convergence
    n_phens = nphns;              // number of phenotypes to check convergence
    epsilon = eps;                // epsilon param to check convergence
    mut_rate = mu;                // mutation rate
    stock_gamete = n_gams;        // number of available gametes
    phen_opt = phenopt;     // optimal phenotype, for for fitness computation
    omega = omg;

//    cout << "Indiv::rdmStart w1 = \n" << w1 << endl;
//    cout << "Indiv::rdmStart w2 = \n" << w2 << endl;
//    cout << "Indiv::rdmStart stock_gamete = \n" << stock_gamete << endl;
}


// Constructor definition
void Indiv::loadParams(mat w_1,
                       mat w_2,
                       int nitr,
                       int nphns,
                       double eps,
                       double mu,
                       int n_gams,
                       rowvec phenopt,
                       double omg)
{
    w1 = w_1;                 // gene-to-gene interaction matrix
    w2 = w_2;                 // gene-to-gene interaction matrix
    niter_conv = nitr;            // number of iterations to check convergence
    n_phens = nphns;              // number of phenotypes to check convergence
    epsilon = eps;                // epsilon param to check convergence
    mut_rate = mu;                // mutation rate
    stock_gamete = n_gams;        // number of available gametes
    phen_opt = phenopt;     // optimal phenotype, for for fitness computation
    omega = omg;                  // omega parameter, for fitness computation
}

// Get networks: To be modified for saving genotypes
mat Indiv::getNetwork1()
{
    return(w1);
}

mat Indiv::getNetwork2()
{
    return(w2);
}

// Get gamete count
int Indiv::getGameteCnt()
{
    return(stock_gamete);
}

// Get phenotype
rowvec Indiv::getPhen()
{
    if( phen_last.n_elem < 1)
    {
        int nrows = w1.n_rows;
        mat A(nrows, nrows);
        A.fill(0.5);
        //cout << "Indiv::getPhen(): w1 = \n" << w1 << endl;
        //cout << "Indiv::getPhen(): A * w1 = \n" << A % w1 << endl;

        // average the two gene networks and load into a network object
        mat w3 = A % w1 + A % w2;

        Network net_ind_tmp;
        net_ind_tmp.loadParams(w3,
        niter_conv,
        n_phens,
        epsilon,
        mut_rate);
        phen_last = net_ind_tmp.getPhen();
        // cout << "Indiv::getPhen(): phen_last = \n" << phen_last << endl;
    }

    return(phen_last);
}

// Get fitness
double Indiv::getFitn()
{
    if( phen_last.n_elem < 1)
    {
        phen_last = getPhen();
    }
    else
    {
        if(phen_last.has_nan())
        {
            fitvalue = 0;
        }
        else
        {
            fitvalue = fitness(phen_last, phen_opt, omega);
        }
    }

    // cout << "getFitn: fitvalue = " << fitvalue << endl;
    // cout << "getFitn: phen_last = " << phen_last << endl;
    // cout << "getFitn: phen_opt = " << phen_opt << endl;
    return(fitvalue);
}

// Produce viable gamete
mat Indiv::getGamete(int infinite)
{
    // cout << "Indiv::getGamete():" << stock_gamete << endl;

    // Recombination events
    // Get matrix dims.
    int n_genes = w1.n_rows;

    // initiate random distributions
    random_device rd;
    mt19937 gen(rd());
    binomial_distribution<> binom_proba(1, 0.5);
    uniform_int_distribution<> unif_gamete(1, 2);

    // Iterate through each gene, and decide if must be swapped by recombination event
    mat w1_rec = w1;
    mat w2_rec = w2;
    for(int loc_idx=0; loc_idx<n_genes; loc_idx++)
    {
        rowvec locus1;
        rowvec locus2;

        bool recomb_event = binom_proba(gen);      // draw random 0 or 1, with 0.5 of success chance
        if(recomb_event > 0)                       // apply recombination if mut_event > 0
        {
            locus1 = w2.row(loc_idx);                              // get locus 1 from w2, and vice-versa
            locus2 = w1.row(loc_idx);                              //
        }
        else
        {
            locus1 = w1.row(loc_idx);                              // get locus 1 from w1, and vice-versa
            locus2 = w2.row(loc_idx);                              //
        }

        w1_rec.row(loc_idx) = locus1;                                     // build recombined matrices
        w2_rec.row(loc_idx) = locus2;
    }

    // return results, decide which gamete must go
    int gamete_trial = unif_gamete(gen);
    if(gamete_trial == 1)
    {
        w_gamete = w1_rec;
    }
    else
    {
        w_gamete = w2_rec;
    }

    // load to Network
    net_gamete.loadParams(w_gamete,
                          niter_conv,
                          n_phens,
                          epsilon,
                          mut_rate);

    // apply mutations
    net_gamete.mutate();

    // check convergence and kill gamete if not convergent
    conv_test = net_gamete.getConvg();


    //cout << "Indiv::getGamete(): phen_gamete = \n" << net_gamete.getPhen() << endl;

    if(conv_test == 0)
    {
        // cout << "Indiv::getGamete() : gamete is not convergent" << endl;
        w_gamete = w_gamete.fill(NAN);
    }
    else
    {
        // cout << "Indiv::getGamete() : gamete is convergent" << endl;
        w_gamete = net_gamete.getNetwork();
    }

    // return output and decrement stock of available gametes, if we are dealing with a standard individual
    if(infinite == 0)
    {
        stock_gamete = stock_gamete - 1;
    }

    return(w_gamete);
}


