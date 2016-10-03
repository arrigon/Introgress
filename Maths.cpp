#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

//rowvec sigmoid(rowvec si)
//{
//    /*
//    ## Sigmoid function, normalizing expression levels between -1 and 1 ##
//    as in Rhone et al 2011, Eq 1, P 2088, sigmoid function f(x) = 2 ⁄ (1 + e^-x)  -1
//     input:
//     - si = raw vector of genes expressions Si
//     output: normalized vector Si
//    */
//
//    rowvec sigm = (2 / (1 + exp(-si))) - 1;
//    return(sigm);
//}

rowvec sigmoid(rowvec si)
{
    /*
    ## Sigmoid function, normalizing expression levels between -1 and 1 ##
    as in Rhone et al 2011, Eq 1, P 2088, sigmoid function f(x) = 2 ⁄ (1 + e^-x)  -1
     input:
     - si = raw vector of genes expressions Si
     output: normalized vector Si
    */

    rowvec sigm = (1 / (1 + exp(-0.4 * si)));
    return(sigm);
}

double sig(mat mat_si, int nlines)
{
    /*
    ## Equilibrium metric, analog to variance measure ##
    Used to assess whether our iterations reached an equilibrium
      as in Siegal & Bergman 2002, P 10530
      inputs:
      - mat_si = matrix of Si values, stored over complete convergence iterations
      - nlines = number of iterations to consider, taken from tail of matrix mat_si
      output: equilibrium metric */
    // correct nlines index
    nlines = nlines - 1;
    int nrows = mat_si.n_rows - 1;

    // get last nlines
    mat mat_si_tmp = mat_si.rows(nrows - nlines, nrows);

    // compute distance to mean
    rowvec mns = mean(mat_si_tmp);
    mat ds = mat_si_tmp.each_row() - mns;
    ds = pow(ds, 2) / (4 * (nlines + 1));
    double ds_sum_mean = as_scalar(mean(sum(ds, 1)));

    // return result
    return(ds_sum_mean);
}

double fitness(rowvec phen, rowvec phen_opt, double omega)
{
    /* Fitness function
    ## as in Rhone et al 2011, Eq 2. F = e^ -((P - Popt)^2 / omega^2)
    */
    double delta = sum(abs(phen - phen_opt));
    double fit = exp(-(pow(delta, 2) / pow(omega, 2)));

    // cout << "fitness: fit = " << fit << "\t delta = " << delta << endl;
    return(fit);
}

rowvec mat2vec(mat w)
{
    // initiate params
    int n_genes = w.n_rows;
    rowvec w_vec(n_genes * n_genes);

    // populate vector
    int cnt = 0;
    for(int i_idx=0; i_idx<n_genes; i_idx++)
    {
        for(int j_idx=0; j_idx<n_genes; j_idx++)
        {
            w_vec.at(cnt) = w.at(i_idx, j_idx);
            cnt = cnt + 1;
        }
    }

    // return output
    return(w_vec);
}

mat vec2mat(rowvec w_vec)
{
    // initiate params
    int n_genes = sqrt(w_vec.n_elem);
    mat w(n_genes, n_genes);

    // populate vector
    int cnt = 0;
    for(int i_idx=0; i_idx<n_genes; i_idx++)
    {
        for(int j_idx=0; j_idx<n_genes; j_idx++)
        {
            w.at(i_idx, j_idx) = w_vec.at(cnt);
            cnt = cnt + 1;
        }
    }

    // return output
    return(w);
}
