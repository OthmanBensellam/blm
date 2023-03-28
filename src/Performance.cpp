#include "Performance.hpp"
#include <cmath>

using namespace std;

Performance::Performance(double T,int degree, int nbTimeSteps, int size, double strike, PnlVect* coeffs) 
{
    T_ = T;
    nbTimeSteps_ = nbTimeSteps;
    size_ = size;
    strike_ = strike;
    degree_ = degree;
    coeffs_ = pnl_vect_create_from_scalar(size,0);
    pnl_vect_clone (coeffs_, coeffs);
}

double Performance::payoff(const PnlMat* path, int tau) {
    double sum = 0;
    double trans;

    for (int d = 0; d < size_; d++) {
        trans = pnl_vect_get(coeffs_, d) * pnl_mat_get(path, d, tau); 
        sum = fmax(sum,trans);
    }
    
    sum -= strike_;

    return fmax(0,sum);
}

double Performance::payoffVect(const PnlVect* St){
    double sum = 0;
    double trans;

    for (int d = 0; d < size_; d++) {
        trans = pnl_vect_get(coeffs_, d) * pnl_vect_get(St, d); 
        sum = fmax(sum,trans);
    }
    
    sum -= strike_;

    return fmax(0,sum);
}