/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This class enforces A and Aww to be trace-free
#ifndef TRACEAREMOVALCARTOON_HPP_
#define TRACEAREMOVALCARTOON_HPP_

//#include "CCZ4Geometry.hpp"
#include "CCZ4CartoonVars.hpp"
#include "Cell.hpp"
#include "Interval.H"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

class TraceARemovalCartoon
{
  public:
    template <class data_t> struct Vars
    {
        Tensor<2, data_t> h;
        Tensor<2, data_t> A;
        data_t hww;
        data_t Aww;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function);
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto vars = current_cell.template load_vars<Vars>();

        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        auto trA = TensorAlgebra::compute_trace(vars.A, h_UU);
        double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);

        trA += (GR_SPACEDIM - CH_SPACEDIM) * vars.Aww / (vars.hww);

        FOR(i, j) { vars.A[i][j] -= one_over_gr_spacedim * vars.h[i][j] * trA; }
        vars.Aww -= one_over_gr_spacedim * vars.hww * trA;

        current_cell.store_vars(vars);
    }
};

template <class data_t>
template <typename mapping_function_t>
void TraceARemovalCartoon::Vars<data_t>::enum_mapping(
    mapping_function_t mapping_function)
{
#if CH_SPACEDIM == 3
    VarsTools::define_symmetric_enum_mapping(mapping_function,
                                             GRInterval<c_h11, c_h33>(), h);
    VarsTools::define_symmetric_enum_mapping(mapping_function,
                                             GRInterval<c_A11, c_A33>(), A);
#elif CH_SPACEDIM == 2
    VarsTools::define_symmetric_enum_mapping(mapping_function,
                                             GRInterval<c_h11, c_h22>(), h);
    VarsTools::define_symmetric_enum_mapping(mapping_function,
                                             GRInterval<c_A11, c_A22>(), A);
#else
#ifdef CH_SPACEDIM
#error define_enum_mapping() has not got your dimension combination implemented.
#endif
#endif
    VarsTools::define_enum_mapping(mapping_function, c_hww, hww);
    VarsTools::define_enum_mapping(mapping_function, c_Aww, Aww);
}

#endif /* TRACEAREMOVALCARTOON_HPP_ */
