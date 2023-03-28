#pragma once

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_array.h"

/// \brief Modèle de Black Scholes
class BlackScholesModel
{
public:
    int size_; /// nombre d'actifs du modèle
    double r_; /// taux d'intérêt
    double rho_; /// paramètre de corrélation
    PnlVect *sigma_; /// vecteur de volatilités
    PnlVect *divid_; /// vecteur des dividendes
    PnlVect *spot_; /// valeurs initiales des sous-jacents
     /**
     * Constructeur de la classe BlackScholes
     */
    BlackScholesModel(int size,double rho,double r, PnlVect* sigma, PnlVect* divid, PnlVect* spot);
    /**
     * Destructeur
     */
    ~BlackScholesModel();
    /**
     * Génère une trajectoire du modèle et la stocke dans path
     *
     * @param[out] path contient une trajectoire du modèle.
     * C'est une matrice de taille (nbTimeSteps+1) x d
     */
    void asset(PnlMat *C, PnlMat *path, double T, double step,int dates, PnlRng *rng);
};
