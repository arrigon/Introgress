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
            // cout << "Network::mutate mutation event" << endl;
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
mat Network::getConvergentNetwork(double netFill, int gradient)
/* Start inititial Network,
   - networkFill [0-1] determines the proportion of non-null interactions in the network
                       N.B. mutations are *not* applied to null interactions, thus preserving the network topology.
   - gradient [0 1, 2 or 3] determines if null interactions are spread randomly in w (gradient = 0)
                       or if one of the following connectivity gradients is applied :
                       gradient = 1 : gradient on rows only
                       gradient = 2 : gradient on columns only
                       gradient = 3 : gradient on rows AND columns
*/
{
    // get params
    networkFillness = netFill;                             // network fillness
    int n_genes = w.n_rows;                                // number of genes
    int n_cells = pow(n_genes, 2);                         // number of cells in w
    int n_nulls = floor(n_cells * (1 - networkFillness));  // number of cells where a zero must be assigned

    // Generate probability gradients of assigning a null value, for each cell of w
    // Start with probas by rows of w
    mat probsRows(n_genes, n_genes);
    for(int gene_idx=0; gene_idx<n_genes; gene_idx++)
    {
        double proba = (gene_idx + 1) /(double) n_genes;   // p = 1/w_row
        rowvec tmp(n_genes);
        tmp.fill(proba);
        probsRows.row(gene_idx) = tmp;
    }

    // set connectivity gradients according to user-defined params
    mat probsCols = probsRows.t();
    mat probsMat;

    // gradient in lines
    if(gradient == 1)
    {
        probsMat = probsRows;
    }

    // gradient in columns
    if(gradient == 2)
    {
        probsMat = probsCols;
    }

    // double gradient : multiply probsRows by ProbsCols
    if(gradient == 3)
    {
        probsMat = probsRows % probsCols;
    }

    // cout << "Network::getConvergentNetwork: connectivity gradient\n" << probsMat << endl;

    // convert to std::vector and feed random sampler with it
    rowvec probsMat_vec = mat2vec(probsMat);          // turn back to std::vector, in order to feed random sampler
    std::vector<double> myProbs;
    for(int idx=0; idx<probsMat_vec.n_elem; idx++)
    {
        double val = probsMat_vec[idx];
        myProbs.push_back(val);
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> rand_int_connect (myProbs.begin(), myProbs.end());


    // Initiate counters
    int conv_ok = 0;
    int myConvgs = 0;

    // Explore random networks until finding a convergent candidate
    while(conv_ok == 0)
    {
        // Initiate random network matrix
        w = randn<mat>(n_genes,n_genes);

        // Assign zeroes in w, on order to decrease the connectivity
        if(gradient > 0)
        {
            // do so according to probability gradients defined earlier
            rowvec zeroes_idx(n_nulls);
            int ok_cnt = 0;
            while(ok_cnt < n_nulls)                    // attribute 0s to w cells, until we reach the amount of null_cells
            {
                int idx = rand_int_connect(gen);       // sample a cell index, using the random sampler initiated earlier (~ proba of null cells)
                uvec check = find(zeroes_idx == idx);  // check if that cell has already been assigned a zero
                if(check.n_elem == 0)                  // if not, proceed further
                {
                    zeroes_idx.at(ok_cnt) = idx;       // save indexes where zero values were assigned
                    w(idx) = 0;                        // assign 0 in w, at focal cell
                    ok_cnt = ok_cnt + 1;               // keep track of how much cells were modified this way
                }
            }
        }
        // cout << "Network::getConvergentNetwork: zeroes_idx " << zeroes_idx << endl;

        if(gradient == 0)
        {
        // assign zeroes randomly in w
        mat w_null = randu<mat>(n_genes,n_genes);
        w.elem( find(w_null > networkFillness)).zeros();
        }



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
    // cout << "Network::getConvergentNetwork: found converging network " << myConvgs << "%" << endl;
    return(w);
}


