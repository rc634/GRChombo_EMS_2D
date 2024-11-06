/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMCONFORMALVARS_HPP_
#define ADMCONFORMALVARS_HPP_

#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

/// Namespace for ADM vars in conformally decomposed form
/** The structs in this namespace collect all the ADM variables. It's main use
 *is to make a local, nicely laid-out, copy of the ADM variables for the
 *current grid cell (Otherwise, this data would only exist on the grid in the
 *huge, flattened Chombo array). \sa {CCZ4Vars, BSSNVars}
 **/
namespace ADMConformalVars
{
/// Vars object for ADM vars, including gauge vars
template <class data_t> struct VarsNoGauge
{
    data_t chi;          //!< Conformal factor
    Tensor<2, data_t> h; //!< Conformal metric
    data_t K;            //!< Trace of the extrinsic curvature
    Tensor<2, data_t> A; //!< trace-free part of the rescale extrinsic
                         //! curvature, i.e. \f$\chi
                         //!(K_{ij})^{\mathrm{TF}}\f$

    // Most sensible to add here rather than BSSN/CCZ4Vars
    data_t phi;
    data_t Pi;
    data_t Ex;
    data_t Ey;
    data_t Ez;
    data_t Bx;
    data_t By;
    data_t Bz;
    data_t Lambda;
    data_t Xi;


    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        // Scalars
        define_enum_mapping(mapping_function, c_chi, chi);
        define_enum_mapping(mapping_function, c_K, K);
        define_enum_mapping(mapping_function, c_phi, phi);
        define_enum_mapping(mapping_function, c_Pi, Pi);
        define_enum_mapping(mapping_function, c_Ex, Ex);
        define_enum_mapping(mapping_function, c_Ey, Ey);
        define_enum_mapping(mapping_function, c_Ez, Ez);
        define_enum_mapping(mapping_function, c_Bx, Bx);
        define_enum_mapping(mapping_function, c_By, By);
        define_enum_mapping(mapping_function, c_Bz, Bz);
        define_enum_mapping(mapping_function, c_Lambda, Lambda);
        define_enum_mapping(mapping_function, c_Xi, Xi);
        // Symmetric 2-tensors
#if CH_SPACEDIM == 3
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_h11, c_h33>(), h);
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_A11, c_A33>(), A);
#elif CH_SPACEDIM == 2
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_h11, c_h22>(), h);
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_A11, c_A22>(), A);
#else
#ifdef CH_SPACEDIM
#error define_enum_mapping() has not got your dimension combination implemented.
#endif
#endif
    }
};

/// Vars object for ADM vars, including gauge vars
template <class data_t> struct VarsWithGauge : public VarsNoGauge<data_t>
{
    data_t lapse;
    Tensor<1, data_t> shift;

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        VarsNoGauge<data_t>::enum_mapping(mapping_function);
        define_enum_mapping(mapping_function, c_lapse, lapse);
#if CH_SPACEDIM == 3
        define_enum_mapping(mapping_function, GRInterval<c_shift1, c_shift3>(),
                            shift);
#elif CH_SPACEDIM == 2
        define_enum_mapping(mapping_function, GRInterval<c_shift1, c_shift2>(),
                            shift);
#else
#ifdef CH_SPACEDIM
#error define_enum_mapping() has not got your dimension combination implemented.
#endif
#endif
    }
};

/// Vars object for ADM vars requiring second derivs, excluding gauge vars
template <class data_t> struct Diff2VarsNoGauge
{
    data_t chi;          //!< Conformal factor
    Tensor<2, data_t> h; //!< Conformal metric
    data_t phi;

    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        define_enum_mapping(mapping_function, c_chi, chi);
        define_enum_mapping(mapping_function, c_phi, phi);
#if CH_SPACEDIM == 3
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_h11, c_h33>(), h);
#elif CH_SPACEDIM == 2
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_h11, c_h22>(), h);
#else
#ifdef CH_SPACEDIM
#error define_enum_mapping() has not got your dimension combination implemented.
#endif
#endif
    }
};

/// Vars object for ADM vars requiring second derivs, with gauge vars
template <class data_t>
struct Diff2VarsWithGauge : public Diff2VarsNoGauge<data_t>
{
    data_t lapse;
    Tensor<1, data_t> shift;

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        Diff2VarsNoGauge<data_t>::enum_mapping(mapping_function);
        define_enum_mapping(mapping_function, c_lapse, lapse);
#if CH_SPACEDIM == 3
        define_enum_mapping(mapping_function, GRInterval<c_shift1, c_shift3>(),
                            shift);
#elif CH_SPACEDIM == 2
        define_enum_mapping(mapping_function, GRInterval<c_shift1, c_shift2>(),
                            shift);
#else
#ifdef CH_SPACEDIM
#error define_enum_mapping() has not got your dimension combination implemented.
#endif
#endif
    }
};
} // namespace ADMConformalVars

#endif /* ADMCONFORMALVARS_HPP_ */
