#include "Basket.hpp"
#include <cmath>
#include <iostream>

using namespace std;

Basket::Basket(double T,int degree, int nbTimeSteps, int size, double strike, PnlVect* coeffs) 
{
    T_ = T;
    nbTimeSteps_ = nbTimeSteps;
    size_ = size;
    strike_ = strike;
    degree_ = degree;
    coeffs_ = pnl_vect_create_from_scalar(size,0);
    pnl_vect_clone (coeffs_, coeffs);
}

double Basket::payoff(const PnlMat* path, int tau) {
    PnlVect* St = pnl_vect_create_from_zero(size_);
    pnl_mat_get_row(St,path,tau);
    double sum = pnl_vect_scalar_prod(coeffs_, St);
    sum -= strike_;

    return fmax(0,sum);
}

double Basket::payoffVect(const PnlVect* St)
{
    double sum = pnl_vect_scalar_prod(coeffs_, St);
    sum -= strike_;

    return fmax(0,sum);
}