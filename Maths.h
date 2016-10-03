// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef MATHS_H
#define MATHS_H

#include <armadillo>
#include <cmath>

arma::rowvec sigmoid(arma::rowvec si);
double sig(arma::mat mat_si, int nlines);
double fitness(arma::rowvec phen, arma::rowvec phen_opt, double omega);
arma::rowvec mat2vec(arma::mat w);
arma::mat vec2mat(arma::rowvec w_vec);
// This is the end of the header guard
#endif
