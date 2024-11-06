/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4VARIABLES_HPP
#define CCZ4VARIABLES_HPP

#include <algorithm>
#include <array>
#include <string>

/// This enum gives the index of the CCZ4 variables on the grid
enum
{
    c_chi,

    c_h11,
    c_h12,
    c_h22,
    c_hww,

    c_K,

    c_A11,
    c_A12,
    c_A22,
    c_Aww,

    c_Theta,

    c_Gamma1,
    c_Gamma2,

    c_lapse,

    c_shift1,
    c_shift2,

    c_B1,
    c_B2,

    c_phi, // scalar field
    c_Pi,                    // scalar field momentum
    c_Lambda,                 // time component of A
    c_Bx,                  // xcomponent of A
    c_By,                  // ycomponent of A
    c_Bz,                  // zcomponent of A
    c_Ex,                  // x electric field
    c_Ey,                  // y electric field
    c_Ez,                  // z electric field
    c_Xi,                  // maxwell constriant
    NUM_CCZ4_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_CCZ4_VARS> ccz4_variable_names = {
    "chi",

    "h11",    "h12",    "h22",

    "hww",

    "K",

    "A11",    "A12",    "A22",

    "Aww",

    "Theta",

    "Gamma1", "Gamma2",

    "lapse",

    "shift1", "shift2",

    "B1",     "B2",

    "phi", "Pi", "Lambda",

    "Bx", "By", "Bz", "Ex", "Ey", "Ez", "Xi"};
} // namespace UserVariables

#endif /* CCZ4VARIABLES_HPP */
