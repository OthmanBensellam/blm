#pragma once

#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"
#include "pnl/pnl_basis.h"

class MonteCarlo
{
public:
    BlackScholesModel *mod_; /// pointeur vers le modèle 
    Option *opt_; /// pointeur sur l'option 
    PnlRng* rng_;/// pointeur sur RNG
    int nbSamples_;/// nombre d'échantillons
    double step_;/// pas de discretisation
    int K_;

     /**
     * Constructeur de la classe MonteCarlo
     */
    MonteCarlo(BlackScholesModel* bs, Option* op, PnlRng* rng, int nbSamples);
    /**
     * Destructeur
     */
    ~MonteCarlo() { }
    /**
     * Calcule le prix de l'option à la date 0
     *
     * @return valeur de l'estimateur Monte Carlo
     */
    void price(double &prix);
    void priceNew(double &prix);
    /**
     * Création de nbSamples path stockés dans l'array de paths Lpath
     * stockage de nbTimeSteps valeurs de St à l'instant tau = T.
     */
    void SetUpPaths(PnlArray* Lpath,  PnlMat* tauPmat, int size, int nbTimeSteps, double T);
    /**
     * Actualisation de matSt (matrice contenant les St)
     * Actualisation de vectPf (vecteur contenant les "payoffs")
     */
    void Actualise(PnlVect* tauIPvect, PnlVect* vectPf,PnlMat* matSt,PnlMat* tauPmat,PnlArray* Lpath,double r, int k);
    /**
     * Choix de l'exercice de l'option
     */
    void Exercise(PnlBasis *B,PnlVect* coef,PnlVect* tauIPvect,PnlMat* matSt,PnlMat* tauPmat, double &meanFinal,double r, int k);

};
