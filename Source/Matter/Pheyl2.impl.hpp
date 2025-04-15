/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #if !defined(PHEYL2_HPP_)
 #error "This file should only be included through Pheyl2.hpp"
 #endif


#ifndef FOR3_LOOP_DEFINED_
#define FOR3_LOOP_DEFINED_
#define FOR1_3D(IDX) for (int IDX = 0; IDX < 3; ++IDX)
#define FOR2_3D(IDX1, IDX2) FOR1_3D(IDX1) FOR1_3D(IDX2)
#define FOR3_3D(IDX1, IDX2, IDX3) FOR2_3D(IDX1, IDX2) FOR1_3D(IDX3)
#define FOR4_3D(IDX1, IDX2, IDX3, IDX4) FOR2_3D(IDX1, IDX2) FOR2_3D(IDX3, IDX4)
#endif


#ifndef PHEYL2_IMPL_HPP_
#define PHEYL2_IMPL_HPP_

template <class data_t> void Pheyl2::compute(Cell<data_t> current_cell) const
{

    // copy data from chombo gridpoint into local variables
    const auto vars = current_cell.template load_vars<Vars>();
    const auto emdvars = current_cell.template load_vars<EMDVars>();

    // Get the coordinates
    const Coordinates<data_t> coords(current_cell, m_dx, m_center);

    // Compute the inverse metric and Christoffel symbols
    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    Tensor<2, data_t, 3> h3_UU = {0.};
    FOR2(i,j) h3_UU[i][j] = h_UU[i][j];
    h3_UU[2][2] = 1./emdvars.hww;

    // Compute the spatial volume element
    const auto epsilon3_LLU = compute_epsilon3_LLU(vars, h3_UU);

    // work out the Newman Penrose scalar
    EMRScalar_t<data_t> out =
        compute_Pheyl2(vars, emdvars, epsilon3_LLU, h3_UU, coords);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Real, c_Pheyl2_Re);
    current_cell.store_vars(out.Im, c_Pheyl2_Im);
}

template <class data_t>
Tensor<3, data_t, 3>
Pheyl2::compute_epsilon3_LLU(const Vars<data_t> &vars,
                            const Tensor<2, data_t, 3> &h3_UU) const
{
    // raised normal vector, NB index 3 is time
    data_t n_U[4];
    n_U[3] = 1. / vars.lapse;
    FOR1(i) { n_U[i] = -vars.shift[i] / vars.lapse; }
    n_U[2] = 0.; // manually set z component (cartoon)

    // 4D levi civita symbol and 3D levi civita tensor in LLL and LUU form
    const auto epsilon4 = TensorAlgebra::epsilon4D();
    Tensor<3, data_t, 3> epsilon3_LLL;
    Tensor<3, data_t, 3> epsilon3_LLU;

    // Projection of antisymmentric Tensor onto hypersurface - see 8.3.17,
    // Alcubierre
    FOR3_3D(i, j, k)
    {
        epsilon3_LLL[i][j][k] = 0.0;
        epsilon3_LLU[i][j][k] = 0.0;
    }
    // projection of 4-antisymetric tensor to 3-tensor on hypersurface
    // note last index contracted as per footnote 86 pg 290 Alcubierre
    FOR3_3D(i, j, k)
    {
        for (int l = 0; l < 4; ++l)
        {
            epsilon3_LLL[i][j][k] += n_U[l] * epsilon4[i][j][k][l] *
                                     vars.lapse / (vars.chi * sqrt(vars.chi));
        }
    }
    // rasing indices
    FOR3_3D(i, j, k)
    {
        FOR1_3D(m)
        {
            epsilon3_LLU[i][j][k] += epsilon3_LLL[i][j][m]
                                       * h3_UU[m][k] * vars.chi;
        }
    }

    return epsilon3_LLU;
}


// Calculation of the Pheyl2 scalar
template <class data_t>
EMRScalar_t<data_t> Pheyl2::compute_Pheyl2(const Vars<data_t> &vars,
                                        const EMDVars<data_t> &emdvars,
                                        const Tensor<3, data_t, 3> &epsilon3_LLU,
                                        const Tensor<2, data_t, 3> &h3_UU,
                                        const Coordinates<data_t> &coords) const
{
    EMRScalar_t<data_t> out;

    // Calculate the tetrads
    const Tetrad_t<data_t> tetrad = compute_null_tetrad(vars, h3_UU, coords);

    // note, u, v, w, upstairs, E, B downstairs
    // u, v, w, are perp to unit normal n_\mu

    // Electromagnetic Fields
    Tensor<1, data_t, 3> E;
    Tensor<1, data_t, 3> B;
    data_t phi = emdvars.phi;
    E[0] = emdvars.Ex;
    E[1] = emdvars.Ey;
    E[2] = emdvars.Ez;
    B[0] = emdvars.Bx;
    B[1] = emdvars.By;
    B[2] = emdvars.Bz;

    // Projection of Electric and Magnetic radiation using tetrads
    out.Real = 0.0;
    out.Im = 0.0;

    // Calculation here
    FOR1_3D(i)
    {
        out.Real += 0.5 * tetrad.v[i] * E[i];
        out.Im += -0.5 * tetrad.w[i] * E[i];
    }
    FOR3_3D(i,j,k)
    {
        out.Real += -0.5 * tetrad.v[i] * tetrad.u[j] * B[k] * epsilon3_LLU[i][j][k];
        out.Im += 0.5 * tetrad.w[i] * tetrad.u[j] * B[k] * epsilon3_LLU[i][j][k];
    }

    data_t f = m_coupling_params.f0
             + m_coupling_params.f1 * phi
             + m_coupling_params.f2 * phi * phi;

    // the -1. in the power is becuase we already square rooted the function
    // out.Real *= exp(-1. * m_coupling_params.alpha * f);
    // out.Im *= exp(-1. * m_coupling_params.alpha * f);

    return out;
}

// Calculation of the null tetrad
// Defintions from gr-qc/0104063
// "The Lazarus project: A pragmatic approach to binary black hole evolutions",
// Baker et al.
template <class data_t>
Tetrad_t<data_t>
Pheyl2::compute_null_tetrad(const Vars<data_t> &vars,
                           const Tensor<2, data_t, 3> &h3_UU,
                           const Coordinates<data_t> &coords) const
{
    Tetrad_t<data_t> out;

    // compute coords
    const data_t x = coords.x;
    const data_t y = coords.y;
    const data_t z = 0.;

    // the alternating levi civita symbol
    const auto epsilon = TensorAlgebra::epsilon_3d_robin();

    // calculate the tetrad
    out.u[0] = x;
    out.u[1] = y;
    out.u[2] = z;

    out.v[0] = -y;
    out.v[1] = x;
    out.v[2] = 0.0;

    out.w[0] = 0.0;
    out.w[1] = 0.0;
    out.w[2] = 0.0;

    // floor on chi
    const data_t chi = simd_max(vars.chi, 1e-4);

    FOR4_3D(i, j, k, m)
    {
        // don't bothr about root chi here, we normalise later
        // we just make w mutually perp to u and v here
        out.w[i] += h3_UU[i][j] * epsilon[j][k][m] * out.v[k] * out.u[m];
    }

    const Tensor<2, data_t, 3> h3 = TensorAlgebra::compute_inverse_sym(h3_UU);

    // Gram Schmitt orthonormalisation
    // Choice of orthonormalisaion to avoid frame-dragging
    data_t omega_11 = 0.0;
    FOR2_3D(i, j) { omega_11 += out.v[i] * out.v[j] * h3[i][j] / chi; }
    FOR1_3D(i) { out.v[i] = out.v[i] / sqrt(omega_11); }

    data_t omega_12 = 0.0;
    FOR2_3D(i, j) { omega_12 += out.v[i] * out.u[j] * h3[i][j] / chi; }
    FOR1_3D(i) { out.u[i] += -omega_12 * out.v[i]; }

    data_t omega_22 = 0.0;
    FOR2_3D(i, j) { omega_22 += out.u[i] * out.u[j] * h3[i][j] / chi; }
    FOR1_3D(i) { out.u[i] = out.u[i] / sqrt(omega_22); }

    data_t omega_13 = 0.0;
    data_t omega_23 = 0.0;
    FOR2_3D(i, j)
    {
        omega_13 += out.v[i] * out.w[j] * h3[i][j] / chi;
        omega_23 += out.u[i] * out.w[j] * h3[i][j] / chi;
    }
    FOR1_3D(i) { out.w[i] += -(omega_13 * out.v[i] + omega_23 * out.u[i]); }

    data_t omega_33 = 0.0;
    FOR2_3D(i, j) { omega_33 += out.w[i] * out.w[j] * h3[i][j] / chi; }
    FOR1_3D(i) { out.w[i] = out.w[i] / sqrt(omega_33); }

    return out;
}

#endif /* PHEYL2_HPP_ */
