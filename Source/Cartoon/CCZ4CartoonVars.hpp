/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4CARTOONVARS_HPP_
#define CCZ4CARTOONVARS_HPP_

#include "ADMConformalVars.hpp"
#include "BSSNVars.hpp"
#include "CCZ4Vars.hpp"
#include "Tensor.hpp"
#include "VarsTools.hpp"

/// Namespace for CCZ4 cartoon vars
/** The structs in this namespace collect all the CCZ4 cartoon variables. It's
 *main use is to make a local, nicely laid-out, copy of the CCZ4 cartoon
 *variables for the current grid cell (Otherwise, this data would only exist on
 *the grid in the huge, flattened Chombo array). \sa {CCZ4Vars,
 *ADMConformalVars}
 **/
namespace CCZ4CartoonVars
{
/// Vars object for CCZ4 cartoon vars, including gauge vars
template <class data_t>
struct VarsNoGauge : public CCZ4Vars::VarsNoGauge<data_t>
{
    data_t hww; //!< Extra-dimensional component of the conformal metric
    data_t Aww; //!< Extra-dimensional component of the traceless part of the
                //!< rescaled extrinsic curvature

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        CCZ4Vars::VarsNoGauge<data_t>::enum_mapping(mapping_function);

        define_enum_mapping(mapping_function, c_hww, hww);
        define_enum_mapping(mapping_function, c_Aww, Aww);
    }
};

/// Vars object for CCZ4 cartoon vars, including gauge vars
template <class data_t>
struct VarsWithGauge : public CCZ4Vars::VarsWithGauge<data_t>
{
    data_t hww; //!< Extra-dimensional component of the conformal metric
    data_t Aww; //!< Extra-dimensional component of the traceless part of the
                //!< rescaled extrinsic curvature

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        CCZ4Vars::VarsWithGauge<data_t>::enum_mapping(mapping_function);
        define_enum_mapping(mapping_function, c_hww, hww);
        define_enum_mapping(mapping_function, c_Aww, Aww);
    }
};

/// Vars object for CCZ4 cartoon vars needing second derivs, excluding gauge
/// vars
template <class data_t>
struct Diff2VarsNoGauge : public ADMConformalVars::Diff2VarsNoGauge<data_t>
{
    data_t hww; //!< Extra-dimensional component of the conformal metric
    data_t Aww; //!< Extra-dimensional component of the traceless part of the
                //!< rescaled extrinsic curvature

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        ADMConformalVars::Diff2VarsNoGauge<data_t>::enum_mapping(
            mapping_function);
        define_enum_mapping(mapping_function, c_hww, hww);
        define_enum_mapping(mapping_function, c_Aww, Aww);
    }
};

/// Vars object for CCZ4 vars needing second derivs, including gauge vars
template <class data_t>
struct Diff2VarsWithGauge : public ADMConformalVars::Diff2VarsWithGauge<data_t>
{
    data_t hww; //!< Extra-dimensional component of the conformal metric
    data_t Aww; //!< Extra-dimensional component of the traceless part of the
                //!< rescaled extrinsic curvature

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        ADMConformalVars::Diff2VarsWithGauge<data_t>::enum_mapping(
            mapping_function);
        define_enum_mapping(mapping_function, c_hww, hww);
        define_enum_mapping(mapping_function, c_Aww, Aww);
    }
};
} // namespace CCZ4CartoonVars

#endif /* CCZ4CARTOONVARS_HPP_ */
