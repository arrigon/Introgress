#include "Maths.h"
#include "Network.h"
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;


// Constructor definition
//Network::Network(mat w_ptr, int nitr, int nphns, double eps, double mu)
//{
//    w = w_ptr;          // gene-to-gene interaction matrix
//    niter_conv = nitr;   // number of iterations to check convergence
//    n_phens = nphns;     // number of phenotypes to check convergence
//    epsilon = eps;       // epsilon param to check convergence
//    mut_rate = mu;       // mutation rate
//    n_genes = w.n_rows;   // Get matrix dims.
//    phen_last = rowvec(n_genes);
//}

void Network::loadParams(mat w_ptr, int nitr, int nphns, double eps, double mu)
{
    w = w_ptr ;          // gene-to-gene interaction matrix
    niter_conv = nitr;   // number of iterations to check convergence
    n_phens = nphns;     // number of phenotypes to check convergence
    epsilon = eps;       // epsilon param to check convergence
    mut_rate = mu;       // mutation rate
    n_genes = w.n_rows;   // Get matrix dims.
    phen_last = rowvec(n_genes);
}

void Network::print()
{
    cout << "Network class parameters:" << endl;
    cout << "w = network of " << w.n_rows << "x" << w.n_rows << " genes" << endl;
    cout << "niter_conv = " << niter_conv << endl;
    cout << "n_phens = " << n_phens << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "mut_rate = " << mut_rate << endl;
}

// Access to convergence test
int Network::getConvg()
{
    convergence();
    return(conv_test);
}

// Access to convergence test
int Network::rerunConvg()
{
    conv_test = -9;
    convergence();
    return(conv_test);
}

// Access to stable phenotype
rowvec Network::getPhen()
{
    convergence();
    return(phen_last);
}

// Access to gene network
mat Network::getNetwork()
{
    return(w);
}

// Apply mutations to w
void Network::mutate()
{
    // initiate random distributions
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> unif_proba(0, 1);
    uniform_int_distribution<> unif_gene(0, n_genes - 1);
    uniform_real_distribution<> unif_effect(-1, 1);

    // Iterate through each gene, and decide if must be modified according to mut_rate
    mat w_mut = w;
    for(int loc_idx=0; loc_idx<n_genes; loc_idx++)
    {
        double mut_event = unif_proba(gen);      // draw random value between 0 and 1
        if(mut_event <= mut_rate)                // apply mutation if mut_event <= mut_rate
        {
            rowvec locus = w.row(loc_idx);                              // extract line
            int mut_pos = unif_gene(gen);                               // get position to mutate
            double mut_size = unif_effect(gen);                         // get mutation effect
            double wij = w_mut(loc_idx, mut_pos);       // modify w matrix accordingly
            if(wij != 0)
            {
                w_mut(loc_idx, mut_pos) = w_mut(loc_idx, mut_pos) + mut_size;       // modify w matrix accordingly
            }
        }
    }

    /* return results
    Note: we are updating the network object directly,
    so we *loose* the original network in the focal object
    */
    w = w_mut;
    conv_test = 0;
    convergence();
}

// Convergence check and computation of phenotype
void Network::convergence()
{
    // initiate variables
    int n_genes = w.n_rows;                                   // get n_genes

    rowvec si = randu<rowvec>(n_genes);       // set initial gene expressions, only positive numbers.
    // si.fill(0.1);
    // cout << "Network::convergence si = \n" << si << endl;

    // iterate through convergence test and store phenotypes along the way
    mat phens = ones<mat>(niter_conv, n_genes);   // set matrix storing phenotype at each iteration
    for(int iter=0; iter<niter_conv; iter++)
    {
        si = si * w;
        si = sigmoid(si);
        phens.row(iter) = si;
    }

    // check convergence and prepare outputs accordingly
    conv_val = sig(phens, n_phens);
    // cout << "Network::convergence conv_val = " << conv_val << endl;

    if(conv_val <= epsilon)
    {
        conv_test = 1;
        phen_last = phens.row(niter_conv - 1);
    }
    else
    {
        conv_test = 0;
        phen_last = phen_last.fill(NAN);
    }

    //cout << "Network::convergence() phen_last =\n" << phen_last << endl;
    //cout << "Network::convergence() conv_val = " << conv_val << endl;
    //cout << "Network::convergence() conv_test = " << conv_test << endl;
}

// Find convergent network
mat Network::getConvergentNetwork(double networkFillness)
/* Start inititial Network, networkFillness [0-1] determines
   the proportion of non-null interactions in the network
   N.B. mutations are *not* applied to null interactions,
   thus preserving the network topology. */
{
    // Initiate counters
    int n_genes = w.n_rows;
    int conv_ok = 0;
    int myConvgs = 0;

    // Explore random networks until finding a convergent candidate
    while(conv_ok == 0)
    {
        // Initiate random network matrix
        w = randn<mat>(n_genes,n_genes);

        // Get null cells
         mat w_null = randu<mat>(n_genes,n_genes);
         w.elem( find(w_null > matrixFillness).zeros();

        // Check convergence of ongoing candidate
        conv_ok = getConvg();

        // If we find one candidate, recheck its stability 100 times
        if(conv_ok == 1)
        {
        myConvgs = 0;
            for(int test_idx=0; test_idx<100; test_idx++)
            {
                int convg2 = rerunConvg();
                myConvgs = myConvgs + convg2;
            }
            if(myConvgs < 100) conv_ok = 0; // discard this candidate if not stable
        }
    }

    // Return outputs
    cout << "main: found converging network " << myConvgs << "%" << endl;
    return(w);
}


