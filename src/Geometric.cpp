#include "Geometric.hpp"
#include <cmath>
#include <iostream>

using namespace std;

Geometric::Geometric(double T, int degree, int nbTimeSteps, int size, double strike) 
{
    T_ = T;
    nbTimeSteps_ = nbTimeSteps;
    size_ = size;
    strike_ = strike;
    degree_ = degree;
}

double Geometric::payoff(const PnlMat* path, int tau) {

    PnlVect* St = pnl_vect_create_from_zero(size_);
    pnl_mat_get_row(St,path,tau);
    double prod = pnl_vect_prod(St);
    
    prod = pow(prod,1/(double)size_);
    prod -= strike_;
    
    return -fmin(0,prod);
}

double Geometric::payoffVect(const PnlVect* St) {

    double prod = pnl_vect_prod(St);
    
    prod = pow(prod,1/(double)size_);
    prod -= strike_;
    
    return -fmin(0,prod);
}