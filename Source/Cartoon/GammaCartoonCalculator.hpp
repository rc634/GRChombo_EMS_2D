/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GAMMACARTOONCALCULATOR_HPP_
#define GAMMACARTOONCALCULATOR_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

class GammaCartoonCalculator
{
    // Only variables needed are metric
    template <class data_t> struct Vars
    {
        Tensor<2, data_t> h;
        data_t hww;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
#if CH_SPACEDIM == 3
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_h11, c_h33>(), h);
#elif CH_SPACEDIM == 2
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_h11, c_h22>(), h);
#else
#ifdef CH_SPACEDIM
#error define_enum_mapping() has not got your dimension combination implemented.
#endif
#endif
            VarsTools::define_enum_mapping(mapping_function, c_hww, hww);
        }
    };

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const double m_dx;
    // std::array<double, CH_SPACEDIM> m_center; //Debugging

  public:
    GammaCartoonCalculator(
        double a_dx /*, std::array<double, CH_SPACEDIM> a_center*/)
        : m_deriv(a_dx), m_dx(a_dx) /*, m_center(a_center)*/
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and calc 1st
        // derivs
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);
        const auto h_UU_ww = 1. / vars.hww;
        const auto chris = compute_christoffel(d1.h, h_UU);
        const int dI = CH_SPACEDIM - 1;
        const int nS =
            GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
        Coordinates<data_t> coords{current_cell, m_dx /*, m_center*/};
#if CH_SPACEDIM == 3
        const double cartoon_coord = coords.z;
#elif CH_SPACEDIM == 2
        const double cartoon_coord = coords.y;
#endif

        const double one_over_cartoon_coord = 1. / cartoon_coord;
        Tensor<1, data_t> chris_ww;
        FOR(i)
        {
            chris_ww[i] = one_over_cartoon_coord *
                          (delta(i, dI) - h_UU[i][dI] * vars.hww);
            FOR(j) chris_ww[i] -= 0.5 * h_UU[i][j] * d1.hww[j];
        }
        Tensor<1, data_t>
            chris_contracted; //!< includes the higher D contributions
        FOR(i)
        {
            chris_contracted[i] =
                chris.contracted[i] + nS * h_UU_ww * chris_ww[i];
        }

#if CH_SPACEDIM == 3
        // assign values of Gamma^k = h_UU^ij * \tilde{Gamma}^k_ij in the output
        // FArrayBox
        current_cell.store_vars(chris_contracted,
                                GRInterval<c_Gamma1, c_Gamma3>());

#elif CH_SPACEDIM == 2
        // assign values of Gamma^k = h_UU^ij * \tilde{Gamma}^k_ij in the output
        // FArrayBox
        current_cell.store_vars(chris_contracted,
                                GRInterval<c_Gamma1, c_Gamma2>());

#else
#ifdef CH_SPACEDIM
#error define_enum_mapping() has not got your dimension combination implemented.
#endif
#endif
    }
};

#endif /* GAMMACARTOONCALCULATOR_HPP_ */
