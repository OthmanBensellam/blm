#pragma once

#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"
#include "pnl/pnl_basis.h"
#include <vector>
#include <functional>
#include </usr/include/python3.10/Python.h>
#include "MonteCarlo.hpp"

class tools
{
public:
    //PnlVect* Coeffs_;
    std::vector<std::function<float(float)>> Base_;
    int m_;
    PnlArray* Lpath ;
    int K_;
    BlackScholesModel* mod_ ;
    MonteCarlo* mc_;
	Option* opt_;        
	int nbSamples_;
	double step_;
    
    tools(MonteCarlo* mc,BlackScholesModel* bs, Option* opt,int nbSamples);
    PnlVect* PolynomeBase(double x);
    double norm(std::function<double(double)> f, int N, int M,PnlMat* Mat);
    double inner_product(std::function<double(double)> f,std::function<double(double)> g, int N, int M, PnlMat* Mat);
    std::vector<std::function<float(float)>> gram_schmidt(int n, int N, int M, PnlMat* Mat);
    int optim();
    double func(double x);
    std::vector<std::function<float(float)>> generateFunctions(int k);
    double integrale(std::function<double(double)> f,double a, double b, int n);
    double Yt(double t );
    double Lambdat (double t,const std::vector<double>& c);
    double Ut(double t,const std::vector<double>& c );
    double Phi(const std::vector<double>& c);
    double TheFct(const std::vector<double>& c);
};
