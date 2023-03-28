#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"

/// \brief Classe Geometric
class Geometric : public Option
{
  public:

    /**
     * Génère un constructeur pour la classe Geometric.
     */
    Geometric(double T, int degree, int nbTimeSteps, int size, double strk);
    /**
     * Destructeur
     */
    ~Geometric() { }
    /**
     * Calcule le payoff de l'option Geometric a l'instant final
     *
     * @return valeur du payoff de l'option Geometric a l'instant final T
     */
    double payoff(const PnlMat* path, int tau);
    /**
     * Calcule le payoff de l'option Geometric pour un vecteur spot St
     *
     * @return valeur du payoff payoff de l'option Geometric pour un vecteur spot St
     */
    double payoffVect(const PnlVect* St);
};
