#include "MonteCarlo.hpp"
#include <cmath>
#include <iostream>
//#include "stlbfgs.h"
#include <cmath>
#include <functional>

using namespace std;

MonteCarlo::MonteCarlo(BlackScholesModel* bs, Option* opt, PnlRng* rng, int nbSamples) {
	mod_ = bs;
	opt_ = opt;
	rng_ = rng;         
	nbSamples_ = nbSamples;
	step_ = opt_->T_/opt_->nbTimeSteps_;
	K_ = 3;
}

void MonteCarlo::SetUpPaths(PnlArray* Lpath,  PnlMat* tauPmat, int size, int nbTimeSteps, double T){

	PnlVect* vectSt = pnl_vect_create_from_scalar(size,0);
	PnlMat* C = pnl_mat_create_from_scalar(size, size, mod_->rho_);

    pnl_mat_set_diag(C, 1, 0);
    pnl_mat_chol(C);

	for (int l = 0; l < nbSamples_; l++){
		PnlMat* path = pnl_mat_create_from_scalar(nbTimeSteps + 1, size, 0); 
        mod_->asset(C,path, T,step_,nbTimeSteps,rng_);
        pnl_mat_get_row(vectSt, path, nbTimeSteps);
        pnl_mat_set_row(tauPmat, vectSt, l);
        pnl_array_set(Lpath, l, (PnlObject*)path);
    }

	pnl_vect_free(&vectSt);
	pnl_mat_free(&C);
}


void MonteCarlo::Actualise(PnlVect* tauIPvect, PnlVect* vectPf,PnlMat* matSt,PnlMat* tauPmat,PnlArray* Lpath, double r, int k){

	PnlVect* vectSt = pnl_vect_create_from_scalar(opt_->size_,0);

	for (int l = 0; l < nbSamples_; l++){
        pnl_mat_get_row (vectSt, tauPmat, l);
		pnl_vect_set (vectPf, l, exp(-r*(pnl_vect_get(tauIPvect,l))*step_)*opt_->payoffVect(vectSt));
		pnl_mat_get_row(vectSt, (PnlMat *)pnl_array_get(Lpath, l), k);
        pnl_mat_set_row(matSt, vectSt, l);
    }

	pnl_vect_free(&vectSt);
}

void MonteCarlo::Exercise(PnlBasis *B,PnlVect* coef,PnlVect* tauIPvect,PnlMat* matSt,PnlMat* tauPmat, double &price, double r, int k){

	PnlVect* vectSt = pnl_vect_create_from_scalar(opt_->size_,0);

	for (int l = 0; l < nbSamples_; l++){
		pnl_mat_get_row (vectSt, matSt, l);
		if(opt_->payoffVect(vectSt)>0){
			if(pnl_basis_eval_vect(B, coef, vectSt)<exp(-r*(k-1)*step_)*opt_->payoffVect(vectSt)){
				pnl_mat_set_row(tauPmat, vectSt, l);
                pnl_vect_set (tauIPvect,l ,k);
			}
        }

		if(k==1){
			pnl_mat_get_row(vectSt,tauPmat, l);
			price += exp(-r*(pnl_vect_get(tauIPvect,l))*step_)*opt_->payoffVect(vectSt);    
        }
    }

	pnl_vect_free(&vectSt);
}

void MonteCarlo::price(double &prix){

	PnlBasis *B = pnl_basis_create_from_degree (PNL_BASIS_HERMITE, opt_->degree_, opt_->size_);
	PnlMat* path = pnl_mat_create_from_scalar(opt_->nbTimeSteps_+1, opt_->size_, 0);
	PnlVect* coef = pnl_vect_create_from_scalar(B->nb_func, 0);
	PnlVect* vectPf = pnl_vect_create_from_scalar(nbSamples_,0);
	PnlMat* matSt = pnl_mat_create_from_scalar(nbSamples_, opt_->size_, 0);
	PnlVect* tauIPvect = pnl_vect_create_from_scalar(nbSamples_,opt_->nbTimeSteps_);
    PnlMat* tauPmat = pnl_mat_create_from_scalar(nbSamples_, opt_->size_, 0);
    PnlArray* Lpath = pnl_array_create(nbSamples_);
	
	int nbTimeSteps = opt_->nbTimeSteps_;
	int size = opt_->size_;
	double r = mod_->r_;
	double T = opt_->T_;
    double alphaH;
	double price =0;

	SetUpPaths(Lpath,  tauPmat, size, nbTimeSteps,T);
	for (int k = nbTimeSteps - 1; k > 0; k--){
		Actualise(tauIPvect, vectPf,matSt,tauPmat,Lpath,r, k);
		pnl_basis_fit_ls (B, coef, matSt, vectPf);
		Exercise(B,coef,tauIPvect,matSt,tauPmat, price,r,k);
	}

	price = price/(double)(nbSamples_);
	prix = max(opt_->payoffVect(mod_->spot_), price);

	pnl_vect_free(&vectPf);
	pnl_vect_free(&coef);
	pnl_vect_free(&tauIPvect);
	pnl_mat_free(&path);
	pnl_mat_free(&matSt);
	pnl_mat_free(&tauPmat);
	pnl_basis_free (&B);
	pnl_array_free(&Lpath);
}





/*
const STLBFGS::func_grad_eval func = [](const std::vector<double> &x, double &f, std::vector<double> &g) {
	f = (x[0] - 7)*(x[0] - 7) +	(x[1] - 1)*(x[1] - 1);
	g[0] = 2*(x[0] - 7);
	g[1] = 2*(x[1] - 1);
};

int MonteCarlo::optim() {
    STLBFGS::Optimizer opt(func);
    std::vector<double> x = {0., 0.};
    opt.run(x);

    std::cout << "Result: x=" << x[0] << ", y=" << x[1] << std::endl;

    if (std::abs(x[0]-7)<1e-3 && std::abs(x[1]-1)<1e-3) {
        std::cout << "Optimization succeeded" << std::endl;
        return 0;
    } else {
        std::cout << "Optimization failed" << std::endl;
        return 1;
    }
}
*/
//use gram shmidt process

/*
int main() {
    const STLBFGS::Optimizer::func_grad_eval func = [](const std::vector<double> &x, double &f, std::vector<double> &g) {
        f = (x[0] - 7)*(x[0] - 7) +
            (x[1] - 1)*(x[1] - 1);
        g[0] = 2*(x[0] - 7);
        g[1] = 2*(x[1] - 1);
    };

    STLBFGS::Optimizer opt{func};
    std::vector<double> x = {10, 10};
    opt.run(x);

    return std::abs(x[0]-7)>1e-3 || std::abs(x[1]-1)>1e-3;
}
*/