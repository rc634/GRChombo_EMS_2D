/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(TWOPUNCTURESINITIALDATA_HPP_)
#error "This file should only be included through TwoPuncturesInitialData.hpp"
#endif

#ifndef TWOPUNCTURESINITIALDATA_IMPL_HPP_
#define TWOPUNCTURESINITIALDATA_IMPL_HPP_

void TwoPuncturesInitialData::compute(Cell<double> current_cell) const
{
    Vars<double> vars;
    // Set only the non-zero components explicitly below
    VarsTools::assign(vars, 0.);

    Coordinates<double> coords(current_cell, m_dx, m_center);
    Tensor<2, double> h_phys, K_tensor;
    Tensor<1, double> shift, Z3;
    double lapse, Theta, hww, Kww;

    interpolate_tp_vars(coords, h_phys, K_tensor, lapse, shift, Theta, Z3, hww,
                        Kww);

    using namespace TensorAlgebra;
    // analytically set Bowen-York properties below (e.g. conformal flatness,
    // tracefree K)

    // metric variables
#if CH_SPACEDIM == 3
    vars.chi = pow(compute_determinant_sym(h_phys), -1. / 3.);
#elif CH_SPACEDIM == 2
    vars.chi = pow(compute_determinant_sym(h_phys) * hww, -1. / 3.);
#endif

    FOR(i, j)
    {
        // Bowen-York data is conformally flat
        vars.h[i][j] = vars.chi * h_phys[i][j];
    }

    // extrinsic curvature
    FOR(i, j) { vars.A[i][j] = vars.chi * K_tensor[i][j]; }
#if CH_SPACEDIM == 3
    // conformal flatness means h_UU = h
    make_trace_free(vars.A, vars.h, vars.h);
#elif CH_SPACEDIM == 2
    vars.hww = vars.chi * hww;
    vars.Aww = vars.chi * Kww;
#endif

    // gauge
    vars.lapse = lapse;

    current_cell.store_vars(vars);
}

void TwoPuncturesInitialData::interpolate_tp_vars(
    const Coordinates<double> &coords, Tensor<2, double> &out_h_phys,
    Tensor<2, double> &out_K_tensor, double &out_lapse,
    Tensor<1, double> &out_shift, double &out_Theta, Tensor<1, double> &out_Z3,
    double &out_hww, double &out_Kww) const
{
    double coords_array[3];
    coords_array[0] = coords.x;
    coords_array[1] = coords.y;

    // If simulation is 3D, fetch values from extra dimension
#if CH_SPACEDIM == 3
    coords_array[2] = coords.z;
#elif CH_SPACEDIM == 2
    coords_array[2] = 0;
#endif

    using namespace TP::Z4VectorShortcuts;
    double TP_state[Qlen];
    m_two_punctures.Interpolate(coords_array, TP_state);

    // metric
    out_h_phys[0][0] = TP_state[g11];
    out_h_phys[0][1] = out_h_phys[1][0] = TP_state[g12];
    out_h_phys[1][1] = TP_state[g22];
#if CH_SPACEDIM == 3
    out_h_phys[0][2] = out_h_phys[2][0] = TP_state[g13];
    out_h_phys[1][2] = out_h_phys[2][1] = TP_state[g23];
    out_h_phys[2][2] = TP_state[g33];
#elif CH_SPACEDIM == 2
    out_hww = TP_state[g33];
#endif

    // extrinsic curvature
    out_K_tensor[0][0] = TP_state[K11];
    out_K_tensor[0][1] = out_K_tensor[1][0] = TP_state[K12];
    out_K_tensor[1][1] = TP_state[K22];
#if CH_SPACEDIM == 3
    out_K_tensor[0][2] = out_K_tensor[2][0] = TP_state[K13];
    out_K_tensor[1][2] = out_K_tensor[2][1] = TP_state[K23];
    out_K_tensor[2][2] = TP_state[K33];
#elif CH_SPACEDIM == 2
    out_Kww = TP_state[K33];
#endif

    // Z4 vector
    out_Z3[0] = TP_state[Z1];
    out_Z3[1] = TP_state[Z2];
    out_Theta = TP_state[Theta];
#if CH_SPACEDIM == 3
    out_Z3[2] = TP_state[Z3];
#endif

    // gauge
    out_lapse = TP_state[lapse];
    out_shift[0] = TP_state[shift1];
    out_shift[1] = TP_state[shift2];
#if CH_SPACEDIM == 3
    out_shift[2] = TP_state[shift3];
#endif
}

#endif /* TWOPUNCTURESINITIALDATA_IMPL_HPP_ */
