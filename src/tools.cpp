#include "tools.hpp"
#include <cmath>
#include <iostream>
//#include "stlbfgs.h"
#include <cmath>
#include <functional>

using namespace std;

//Implementation call simple (un seul sous-jacent)

/*
PnlVect* tools::PolynomeBase(double x){
	PnlVect* PolyVect = pnl_vect_create_from_scalar(K_+1,0.);
	for(int i=0;i<K_+1;i++){
		pnl_vect_set(PolyVect,i,pow(x,i));
	}
	return PolyVect;
}
*/

tools::tools(MonteCarlo* mc, BlackScholesModel* bs, Option* opt,int nbSamples) {
    
    //std::vector<std::function<float(float)>> Base_;
    PnlArray* Lpath = pnl_array_create(nbSamples);
    PnlMat* tauPmat = pnl_mat_create_from_scalar(nbSamples_, opt_->size_, 0);
    mc -> SetUpPaths(Lpath,  tauPmat, opt_->size_, opt_->nbTimeSteps_, opt_->T_);
    mc_ = mc;
	mod_ = bs;
	opt_ = opt;        
	nbSamples_ = nbSamples;
	step_ = opt_->T_/opt_->nbTimeSteps_;
	K_ = 5;
    m_ = 0;
    //PnlVect* Coeffs_ = pnl_vect_create_from_scalar(K_,1.0);
}

double tools::norm(std::function<double(double)> f, int N, int M, PnlMat* Mat){
	double result;
	for(int n =0;n<N;n++){
		for(int m = 0; m<M;m++){
			result += f(pnl_mat_get(Mat,n,m));
		}
	}
	return result/(double)(N*M);
}

double tools::inner_product(std::function<double(double)> f,std::function<double(double)> g, int N, int M, PnlMat* Mat) {
   // Define the integrand for computing the inner product
   auto integrand = [f,g](double x) { return f(x) * g(x); };
   
   // Use the polarization identity to compute the inner product
   double u_plus_v_norm_sq = norm([f,g](double x) { return f(x) + g(x); }, N, M, Mat);
   double u_minus_v_norm_sq = norm([f,g](double x) { return f(x) - g(x); }, N, M, Mat);
   double result = (u_plus_v_norm_sq * u_plus_v_norm_sq - u_minus_v_norm_sq * u_minus_v_norm_sq) / 4.0;
   
   return result;
}

std::vector<std::function<float(float)>> tools::generateFunctions(int k) {
    std::vector<std::function<float(float)>> functions;
    for (int i = 0; i <= k; i++) {
        functions.push_back([i](float x) {
            float result = 1;
            for (int j = 0; j < i; j++) {
                result *= x;
            }
            return result;
        });
    }
    return functions;
}


std::vector<std::function<float(float)>> tools::gram_schmidt(int n, int N, int M, PnlMat* Mat) {
	std::vector<std::function<float(float)>> Base = generateFunctions(K_);
    for (int i = 0; i < n; i++) {
        // Subtract off the projection onto previous basis vectors
        for (int j = 0; j < i; j++) {
            double proj = inner_product(Base[j], Base[i],N, M, Mat) / inner_product(Base[i], Base[i],N, M, Mat);
			Base[i] = ([i,j,Base,proj](float x){return Base[i](x)-proj*Base[j](x);});
            //basis[i] -= proj*basis[j];
        }
        // Normalize the basis vector
        double nrm = sqrt(inner_product(Base[i], Base[i],N, M, Mat));
		Base[i] = ([i,Base,nrm](float x){return Base[i](x)/(double)nrm;});
		//pnl_vect_set(basis,i,pnl_vect_get(basis,i)/(double)nrm);
    }
	return Base;
}

//remplacer par la fonction native c++ 
double tools::integrale(std::function<double(double)> f,double a, double b, int n) {
    // Calcule l'intégrale de f(x) de a à b en utilisant l'interpolation constante par morceaux avec n subdivisions
    double dx = (b - a) / n;
    double x = a;
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double mid = (x + x + dx) / 2.0;
        sum += f(mid) * dx;
        x += dx;
    }
    return sum;
}

double tools::Yt(double t ){
    PnlVect* vectSt = pnl_vect_create_from_scalar(opt_->size_,0);
    int t_index = (int) t/step_;
    pnl_mat_get_row(vectSt, (PnlMat *)pnl_array_get(Lpath, m_), t_index);
    return exp(-t*mod_->r_)*opt_->payoffVect(vectSt);
}

//doit etre modifié pr atetidre le multi-jacents.
double tools::Lambdat(double t,const std::vector<double>& c ){
    PnlVect* vectSt = pnl_vect_create_from_scalar(opt_->size_,0);
    int t_index = (int) t/step_;
    pnl_mat_get_row(vectSt, (PnlMat *)pnl_array_get(Lpath, m_), t_index);
    if(opt_->payoffVect(vectSt)<=0){
        return 0;
    }
    //calcule interieur du polynome
    PnlVect* vectStLOG = pnl_vect_create_from_scalar(opt_->size_,0);
    //pnl_vect_map (vectStLOG, vectSt, log);
    double p = 0;
    for(int j =0; j<opt_->size_; j++){
        for(int i =0; i<K_+1 ; i++){
            p += c[i]*Base_[i](pnl_vect_get(vectStLOG,j));
        }
    }
    return exp(p);
}

double tools::Ut(double t,const std::vector<double>& c){
    //calcule de l'integrale
    int n = 1000;//choisir n conveable 
    double dx = (t) / (double)n;
    double x = 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double mid = (x + x + dx) / 2.0;
        sum += Lambdat(mid,c) * dx;
        x += dx;
    }
    return exp(-sum);
}

double tools::Phi(const std::vector<double>& c){
    //Ut*Yt
    double result = Ut(opt_->T_,c)*Yt(opt_->T_);
    //l'exp intgrl
    int n = 1000;//choisir n conveable 
    double dx = (opt_->T_) / (double)n;
    double x = 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double mid = (x + x + dx) / 2.0;
        sum += Ut(mid,c)*Yt(mid)*Lambdat(mid,c) * dx;
        x += dx;
    }
    return result + sum;
}

double tools::TheFct(const std::vector<double>& c){
    double result = 0;
    while(m_< nbSamples_){
        result += Phi(c);
        m_++;
    }
    return result/(double)nbSamples_;
}

int tools::optim(){
        // Initialisation de l'interpréteur Python
    Py_Initialize();

    // Importation du module scipy.optimize
    PyObject* pModule = PyImport_ImportModule("scipy.optimize");

    // Récupération de la fonction minimize
    PyObject* pFunc = PyObject_GetAttrString(pModule, "minimize");

    // Création des arguments de la fonction minimize
    PyObject* pArgs = PyTuple_New(2);

    // Création de la fonction à minimiser
    PyObject* pFuncCpp = PyCapsule_New((void*)&tools::TheFct, "my_function", NULL);

    // Création du vecteur d'initialisation
    PyObject* pInit = PyList_New(K_+1);
    for(int i =0; i<K_+1; i++){
        PyList_SET_ITEM(pInit, i, PyFloat_FromDouble(1.0));
    }

    // Configuration des arguments de la fonction minimize
    PyTuple_SET_ITEM(pArgs, 0, pFuncCpp);
    PyTuple_SET_ITEM(pArgs, 1, pInit);

    // Appel de la fonction minimize
    PyObject* pResult = PyObject_CallObject(pFunc, pArgs);

    // Récupération de la valeur minimale
    double minimum = PyFloat_AsDouble(PyObject_GetAttrString(pResult, "fun"));

    // Libération de la mémoire
    Py_DECREF(pModule);
    Py_DECREF(pFunc);
    Py_DECREF(pArgs);
    Py_DECREF(pFuncCpp);
    Py_DECREF(pInit);
    Py_DECREF(pResult);

    // Fermeture de l'interpréteur Python
    Py_Finalize();

    std::cout << "Minimum : " << minimum << std::endl;

    return 0;
}


//C'est la qu'on maximise par rapports qux Coeffs_



/*
const STLBFGS::func_grad_eval func = [](const std::vector<double> &x, double &f, std::vector<double> &g) {
	f = (x[0] - 7)*(x[0] - 7) +	(x[1] - 1)*(x[1] - 1);
	g[0] = 2*(x[0] - 7);
	g[1] = 2*(x[1] - 1);
};

int tools::optim() {
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
