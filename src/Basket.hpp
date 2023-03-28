#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"

/// \brief Classe Basket
class Basket : public Option
{
  public:
    PnlVect* coeffs_;///Coefficients de ponderation des sous-jacents.

    /**
     *Génère un constructeur pour la classe Basket
     */
    Basket(double T, int degree, int nbTimeSteps, int size, double strk,PnlVect* coeffs);
    /**
     * Destructeur
     */
    ~Basket() { }
    /**
     * Calcule le payoff de l'option Basket a l'instant final
     *
     * @return valeur du payoff de l'option Basket a l'instant final T
     */
    double payoff(const PnlMat* path, int tau);
    /**
     * Calcule le payoff de l'option Basket pour un vecteur spot St
     *
     * @return valeur du payoff payoff de l'option Basket pour un vecteur spot St
     */
    double payoffVect(const PnlVect* St);
};