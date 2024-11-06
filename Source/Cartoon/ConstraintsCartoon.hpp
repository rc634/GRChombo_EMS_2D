/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This compute class calculates Hamiltonian and Momentum constraints for
// cartoon reduced equations

#ifndef CONSTRAINTSCARTOON_HPP_
#define CONSTRAINTSCARTOON_HPP_

#include "CCZ4CartoonVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FArrayBox.H"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "simd.hpp"

#include "CCZ4CartoonGeometry.hpp"

#include <array>

template <class potential_t> class Constraints
{
  public:
    /// CCZ4 variables
    template <class data_t> using Vars = CCZ4CartoonVars::VarsNoGauge<data_t>;

    /// CCZ4 variables
    template <class data_t>
    using Diff2Vars = CCZ4CartoonVars::Diff2VarsNoGauge<data_t>;

    template <class data_t> struct constraints_t
    {
        data_t Ham;
        Tensor<1, data_t> Mom;
        data_t Mom_abs;

        data_t rho;
        data_t rho_ADM;
        Tensor<1, data_t> Si;
        Tensor<2, data_t> Sij;
        data_t Sww;
        data_t S;
        Tensor<2, data_t> Sij_TF;
        data_t Sww_TF;
    };

    Constraints(double dx, potential_t a_potential, double a_G_Newton,
                double cosmological_constant = 0);

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const FourthOrderDerivatives m_deriv;
    double m_cosmological_constant;
    double m_dx;
    potential_t m_potential;
    double m_G_Newton;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    constraints_t<data_t>
    constraint_equations(const vars_t<data_t> &vars,
                         const vars_t<Tensor<1, data_t>> &d1,
                         const diff2_vars_t<Tensor<2, data_t>> &d2,
                         const double &cartoon_coord) const;
};

#include "ConstraintsCartoon.impl.hpp"

#endif /* CONSTRAINTSCARTOON_HPP_ */
