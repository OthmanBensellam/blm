#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Classe Option abstraite
class Option
{
public:
    double T_; /// maturité
    double strike_;///Strike
    int dates_; /// nombre de dates d'exercice
    int size_; /// dimension du modèle, redondant avec BlackScholesModel::size_
    int nbTimeSteps_;///nombre de pas de descritisation
    int degree_;///degree
    /**
     * Calcule le payoff de l'option Basket a l'instant final
     */
    virtual double payoff(const PnlMat *path, int tau) = 0;
    /**
     * Calcule le payoff de l'option Basket pour a vecteur spot St
     */
    virtual double payoffVect(const PnlVect *St) = 0;
};


