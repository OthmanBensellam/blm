#include "BlackScholesModel.hpp"
#include <cmath>
#include <iostream>

using namespace std;

BlackScholesModel::BlackScholesModel(int size,double rho,double r, PnlVect* sigma, PnlVect* divid, PnlVect* spot){

    size_ = size;
    rho_ = rho;
    r_ = r;
    sigma_ = pnl_vect_create_from_scalar(size,0);
    pnl_vect_clone (sigma_, sigma);
    divid_ = pnl_vect_create_from_scalar(size,0);
    pnl_vect_clone (divid_, divid);
    spot_ = pnl_vect_create_from_scalar(size,0);
    pnl_vect_clone (spot_, spot);
}

BlackScholesModel::~BlackScholesModel()
{
    pnl_vect_free(&sigma_);
	pnl_vect_free(&divid_);
    pnl_vect_free(&spot_);
}

void BlackScholesModel::asset(PnlMat * C, PnlMat *path, double T, double step,int nbTimeSteps, PnlRng *rng){

    PnlVect* Gk = pnl_vect_create(size_);
	PnlVect* L = pnl_vect_create(size_);
    double div;
	double sigma;
    double expon;
    double value;

    pnl_mat_set_row(path, spot_, 0);
	for (int k = 1; k < nbTimeSteps + 1; k++){
		pnl_vect_rng_normal(Gk, size_, rng);
		for (int d = 0; d < size_; d++){
            div = pnl_vect_get(divid_,d);
            sigma = pnl_vect_get(sigma_, d);
			pnl_mat_get_row(L, C, d);
            expon = exp((r_ - div - sigma * sigma / (double)2) * step + sigma * sqrt(step) * pnl_vect_scalar_prod(L, Gk));
            value = expon * pnl_mat_get(path,k-1,d);
            pnl_mat_set(path,k,d,value);
		}
	}
    pnl_vect_free(&L);
    pnl_vect_free(&Gk);

}