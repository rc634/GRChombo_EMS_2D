/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXPERIMENTALGAUGE_HPP_
#define EXPERIMENTALGAUGE_HPP_

#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "MovingPunctureGauge.hpp"
#include "Tensor.hpp"

/// This is an example of a gauge class that can be used in the CCZ4RHS compute
/// class
/**
 * robins experimental gauge to try fix EMS boosts
 **/
class ExperimentalGauge
{
  public:
    using params_t = MovingPunctureGauge::params_t;

    // the decay parameter
    double m_kappa = 1.;

  protected:
    params_t m_params;

    /// Vars needed internally in 'compute'
    template <class data_t> struct Vars
    {
        Tensor<1, data_t> shift;
        Tensor<1, data_t> Gamma; //!< Conformal connection functions
        data_t chi;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_shift1, c_shift2>(), shift);
            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_Gamma1, c_Gamma2>(), Gamma);
            VarsTools::define_enum_mapping( mapping_function, c_chi, chi);
        }
    };

  public:
    ExperimentalGauge(const params_t &a_params) : m_params(a_params)
    {
    }

    // compute function template
    // can be used for initial data or post time step stuff ..
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();

        Tensor<1, data_t> B;
        Tensor<1, data_t> shift;
        FOR(i)
        {
            //shift[i] = vars.chi * vars.shift[i];
            shift[i] = vars.shift[i];
            B[i] = m_params.shift_Gamma_coeff * vars.Gamma[i] -
                   m_params.eta * shift[i];
        }

        current_cell.store_vars(B, GRInterval<c_B1, c_B2>());
        current_cell.store_vars(shift, GRInterval<c_shift1, c_shift2>());
    }

    // the gauge evolution equations
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    inline void rhs_gauge(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                          const vars_t<Tensor<1, data_t>> &d1,
                          const diff2_vars_t<Tensor<2, data_t>> &d2,
                          const vars_t<data_t> &advec) const
    {
        rhs.lapse = 1. * advec.lapse
                     - 2. * vars.lapse * (vars.K - 2 * vars.Theta); // there was a 1.7 coeff here
        FOR(i)
        {
            rhs.shift[i] = 1. * advec.shift[i] +
                           m_params.shift_Gamma_coeff * vars.Gamma[i] -
                           m_params.eta * vars.shift[i] - vars.B[i];
            rhs.B[i] = - 0.1 * vars.B[i] ; // static, stays the same to save initial condition
        }
    }
};

#endif /* EXPERIMENTALGAUGE_HPP_ */
