/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BLACKSTRINGTAGGINGCRITERION_HPP_
#define BLACKSTRINGTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <math.h>

class BlackStringTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_chiTagMultiplier;
    const double m_KTagMultiplier;
    const double m_gammaTagMultiplier;
    const FourthOrderDerivatives m_deriv;

  public:
    BlackStringTaggingCriterion(double dx, double chiTagMultiplier,
                                double KTagMultiplier,
                                double gammaTagMultiplier)
        : m_dx(dx), m_chiTagMultiplier(chiTagMultiplier),
          m_KTagMultiplier(KTagMultiplier),
          m_gammaTagMultiplier(gammaTagMultiplier), m_deriv(dx){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto chi = current_cell.load_vars(c_chi);
        Tensor<1, data_t> d1_chi;
        FOR(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        data_t mod_d1_chi = 0;
        FOR(idir) mod_d1_chi += d1_chi[idir] * d1_chi[idir];

        data_t criterion =
            m_chiTagMultiplier * m_dx * sqrt(mod_d1_chi) / pow(chi, 2);

        // Add tagging based on K for the beginning of the simulation
        auto K = current_cell.load_vars(c_K);
        criterion += m_KTagMultiplier * m_dx * abs(K);

        Tensor<1, data_t> Gamma;
        Gamma[0] = current_cell.load_vars(c_Gamma1);
        Gamma[1] = current_cell.load_vars(c_Gamma2);
        criterion += m_gammaTagMultiplier * m_dx *
                     sqrt(pow(Gamma[0], 2) + pow(Gamma[1], 2));

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* BLACKSTRINGTAGGINGCRITERION_HPP_ */
