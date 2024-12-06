/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(WEYL4_HPP_)
#error "This file should only be included through Weyl4.hpp"
#endif

// Initialize new 3D variables to pass into EB function.
#ifndef FOR3_LOOP_DEFINED_WEYL_
#define FOR3_LOOP_DEFINED_WEYL_
#define FOR_3D(IDX) for (int IDX = 0; IDX < 3; ++IDX)
#define FOR_3D(IDX1, IDX2) FOR_3D(IDX1) FOR_3D(IDX2)
#define FOR_3D(IDX1, IDX2, IDX3) FOR_3D(IDX1, IDX2) FOR_3D(IDX3)
#define FOR_3D(IDX1, IDX2, IDX3, IDX4) FOR_3D(IDX1, IDX2) FOR_3D(IDX3, IDX4)
#endif


#ifndef WEYL4_IMPL_HPP_
#define WEYL4_IMPL_HPP_

// ******** 2D code Specific ********
// Since we have only h11, h12, h22 and hww in our problem,
// we need to unfold the components to get the full 3x3 matrices
// The Stragtegy is to load the CCZ4Cartoon variables in,
// reconstruct the normal CCZ4Cartoon variables in 3D,
// then use those to compute EB fields and Psi4.



template <class data_t>
CCZ4vars3D<data_t> Weyl4::load_ccz4(const Vars<data_t> &vars,
                                    const Vars<Tensor<1, data_t>> &d1,
                                    const Diff2Vars<Tensor<2, data_t>> &d2,
                                    const Coordinates<data_t> &coords) const
{
    // The z here is radial axis of a normal cynlindrical coordinate.
    const double z = coords.y;

    CCZ4vars3D<data_t> output;

    output.lapse = 0.0;
    output.chi = 0.0;
    output.K = 0.0;

    FOR_3D(i)
    {
        output.shift[i] = 0.0;
        output.Gamma[i] = 0.0;
        output.B[i] = 0.0;
        output.d1_lapse[i] = 0.0;
        output.d1_chi[i] = 0.0;
        output.d1_K[i] = 0.0;
    }

    FOR_3D(i, j)
    {
        output.h[i][j] = 0.0;
        output.A[i][j] = 0.0;
        output.d1_shift[i][j] = 0.0;
        output.d1_Gamma[i][j] = 0.0;
        output.d1_B[i][j] = 0.0;
        output.d2_chi[i][j] = 0.0;
    }

    FOR_3D(i, j, k)
    {
        output.d1_h[i][j][k] = 0.0;
        output.d1_A[i][j][k] = 0.0;
    }

    FOR_3D(i, j, k, l) { output.d2_h[i][j][k][l] = 0.0; }

    // Copy trivials component

    output.lapse = vars.lapse;
    output.chi = vars.chi;
    output.K = vars.K;

    FOR(i)
    {
        output.shift[i] = vars.shift[i];
        output.Gamma[i] = vars.Gamma[i];
        output.B[i] = vars.B[i];
        output.d1_lapse[i] = d1.lapse[i];
        output.d1_chi[i] = d1.chi[i];
        output.d1_K[i] = d1.K[i];
    }

    FOR(i, j)
    {
        output.h[i][j] = vars.h[i][j];
        output.A[i][j] = vars.A[i][j];
        output.d1_shift[i][j] = d1.shift[i][j];
        output.d1_Gamma[i][j] = d1.Gamma[i][j];
        output.d1_B[i][j] = d1.B[i][j];
        output.d2_chi[i][j] = d2.chi[i][j];
    }
    FOR(i, j, k)
    {
        output.d1_h[i][j][k] = d1.h[i][j][k];
        output.d1_A[i][j][k] = d1.A[i][j][k];
    }

    FOR(i, j, k, l) { output.d2_h[i][j][k][l] = d2.h[i][j][k][l]; }

    // Assuming lapse,chi,K are scalar densities while shift,Gamma and B are
    // vector

    // Filling vars

    output.h[2][2] = vars.hww;
    output.A[2][2] = vars.Aww;

    // Filling d1

    // First derivatives of related scalar densities are all 0

    // Vector densities
    output.d1_shift[2][2] = vars.shift[1] / z;
    output.d1_Gamma[2][2] = vars.Gamma[1] / z;
    output.d1_B[2][2] = vars.B[1] / z;

    // tensor densities
    output.d1_h[0][2][2] = vars.h[0][1] / z;
    output.d1_h[2][0][2] = vars.h[1][0] / z;
    output.d1_h[1][2][2] = (vars.h[1][1] - vars.hww) / z;
    output.d1_h[2][1][2] = (vars.h[1][1] - vars.hww) / z;
    output.d1_A[0][2][2] = vars.A[0][1] / z;
    output.d1_A[2][0][2] = vars.A[1][0] / z;
    output.d1_A[1][2][2] = (vars.A[1][1] - vars.Aww) / z;
    output.d1_A[2][1][2] = (vars.A[1][1] - vars.Aww) / z;

    // Filling d2

    // Scalar densities

    output.d2_chi[2][2] = d1.chi[1] / z;

    // Vector densities

    // Tensor densities

    output.d2_h[2][2][2][2] =
        2 * (vars.h[1][1] - vars.hww) / z / z + d1.hww[1] / z;

    output.d2_h[0][2][0][2] = d1.h[0][1][0] / z;
    output.d2_h[1][2][0][2] = (d1.h[1][1][0] - d1.hww[0]) / z;
    output.d2_h[0][2][1][2] = d1.h[0][1][1] / z - vars.h[0][1] / z / z;
    output.d2_h[1][2][1][2] =
        (d1.h[1][1][1] - d1.hww[1]) / z - (vars.h[1][1] - vars.hww) / z / z;

    output.d2_h[0][2][2][0] = d1.h[0][1][0] / z;
    output.d2_h[1][2][2][0] = (d1.h[1][1][0] - d1.hww[0]) / z;
    output.d2_h[0][2][2][1] = d1.h[0][1][1] / z - vars.h[0][1] / z / z;
    output.d2_h[1][2][2][1] =
        (d1.h[1][1][1] - d1.hww[1]) / z - (vars.h[1][1] - vars.hww) / z / z;

    output.d2_h[2][0][0][2] = d1.h[0][1][0] / z;
    output.d2_h[2][1][0][2] = (d1.h[1][1][0] - d1.hww[0]) / z;
    output.d2_h[2][0][1][2] = d1.h[0][1][1] / z - vars.h[0][1] / z / z;
    output.d2_h[2][1][1][2] =
        (d1.h[1][1][1] - d1.hww[1]) / z - (vars.h[1][1] - vars.hww) / z / z;

    output.d2_h[2][0][2][0] = d1.h[0][1][0] / z;
    output.d2_h[2][1][2][0] = (d1.h[1][1][0] - d1.hww[0]) / z;
    output.d2_h[2][0][2][1] = d1.h[0][1][1] / z - vars.h[0][1] / z / z;
    output.d2_h[2][1][2][1] =
        (d1.h[1][1][1] - d1.hww[1]) / z - (vars.h[1][1] - vars.hww) / z / z;

    output.d2_h[0][0][2][2] = d1.h[0][0][1] / z;
    output.d2_h[0][1][2][2] = d1.h[0][1][1] / z - (vars.h[1][0]) / z / z;
    output.d2_h[1][0][2][2] = d1.h[1][0][1] / z - (vars.h[1][0]) / z / z;
    output.d2_h[1][1][2][2] =
        d1.h[1][1][1] / z - 2 * (vars.h[1][1] - vars.hww) / z / z;

    return output;

    // ******** End of 2D code specific transforms ********
}

template <class data_t> void Weyl4::compute(Cell<data_t> current_cell) const
{

    // copy data from chombo gridpoint into local variables
    auto vars = current_cell.template load_vars<Vars>();
    auto d1 = m_deriv.template diff1<Vars>(current_cell);
    auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

    // Get the coordinates
    const Coordinates<data_t> coords(current_cell, m_dx, m_center);

    // Transform Cartoon variables to traditional 3D Chombo variables

    CCZ4vars3D<data_t> vars_ccz43D = load_ccz4(vars, d1, d2, coords);

    // Compute the E and B fields
    EBFields_t<data_t> ebfields = compute_EB_fields(vars_ccz43D, coords);

    // work out the Newman Penrose scalar
    NPScalar_t<data_t> out = compute_Weyl4(ebfields, vars_ccz43D, coords);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Real, c_Weyl4_Re);
    current_cell.store_vars(out.Im, c_Weyl4_Im);
}

// Calculation of E and B fields, using tetrads from gr-qc/0104063
// Formalism from Alcubierre book
template <class data_t>
EBFields_t<data_t>
Weyl4::compute_EB_fields(const CCZ4vars3D<data_t> &vars,
                         const Coordinates<data_t> &coords) const
{
    EBFields_t<data_t> out;

    // raised normal vector, NB index 3 is time
    data_t n_U[4];
    n_U[3] = 1. / vars.lapse;
    FOR_3D(i) { n_U[i] = -vars.shift[i] / vars.lapse; }

    // 4D levi civita symbol and 3D levi civita tensor in LLL and LUU form
    const std::array<std::array<std::array<std::array<int, 4>, 4>, 4>, 4>
        epsilon4 = TensorAlgebra::epsilon4D();
    Tensor<3, data_t, 3> epsilon3_LLL;
    Tensor<3, data_t, 3> epsilon3_LUU;

    // Projection of antisymmentric Tensor onto hypersurface - see 8.3.17,
    // Alcubierre
    FOR_3D(i, j, k)
    {
        epsilon3_LLL[i][j][k] = 0.0;
        epsilon3_LUU[i][j][k] = 0.0;
    }
    // projection of 4-antisymetric tensor to 3-tensor on hypersurface
    // note last index contracted as per footnote 86 pg 290 Alcubierre
    FOR_3D(i, j, k)
    {
        for (int l = 0; l < 4; ++l)
        {
            epsilon3_LLL[i][j][k] += n_U[l] * epsilon4[i][j][k][l] *
                                     vars.lapse * pow(vars.chi, -1.5);
        }
    }
    // rasing indices
    data_t deth = vars.h[0][0] * vars.h[1][1] * vars.h[2][2] +
                  2 * vars.h[0][1] * vars.h[0][2] * vars.h[1][2] -
                  vars.h[0][0] * vars.h[1][2] * vars.h[1][2] -
                  vars.h[1][1] * vars.h[0][2] * vars.h[0][2] -
                  vars.h[2][2] * vars.h[0][1] * vars.h[0][1];

    data_t deth_inverse = 1. / deth;
    Tensor<2, data_t, 3> h_UU;
    h_UU[0][0] = (vars.h[1][1] * vars.h[2][2] - vars.h[1][2] * vars.h[1][2]) *
                 deth_inverse;
    h_UU[0][1] = (vars.h[0][2] * vars.h[1][2] - vars.h[0][1] * vars.h[2][2]) *
                 deth_inverse;
    h_UU[0][2] = (vars.h[0][1] * vars.h[1][2] - vars.h[0][2] * vars.h[1][1]) *
                 deth_inverse;
    h_UU[1][1] = (vars.h[0][0] * vars.h[2][2] - vars.h[0][2] * vars.h[0][2]) *
                 deth_inverse;
    h_UU[1][2] = (vars.h[0][1] * vars.h[0][2] - vars.h[0][0] * vars.h[1][2]) *
                 deth_inverse;
    h_UU[2][2] = (vars.h[0][0] * vars.h[1][1] - vars.h[0][1] * vars.h[0][1]) *
                 deth_inverse;
    h_UU[1][0] = h_UU[0][1];
    h_UU[2][0] = h_UU[0][2];
    h_UU[2][1] = h_UU[1][2];

    FOR_3D(i, j, k)
    {
        FOR_3D(m, n)
        {
            epsilon3_LUU[i][j][k] += epsilon3_LLL[i][m][n] * h_UU[m][j] *
                                     vars.chi * h_UU[n][k] * vars.chi;
        }
    }

    // Extrinsic curvature
    Tensor<2, data_t, 3> K_tensor;
    Tensor<3, data_t, 3> d1_K_tensor;
    Tensor<3, data_t, 3> covariant_deriv_K_tensor;

    // Compute inverse, Christoffel symbols and Ricci Tensor
    using namespace TensorAlgebra;
    //    const auto chris = compute_christoffel(vars.d1_h, h_UU);
    //    const auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU,
    //    chris);

    // Compute full spatial Christoffel symbols
    // const Tensor<3, data_t,3> chris_phys =
    //     compute_phys_chris(vars.d1_chi, vars.chi, vars.h, h_UU, chris.ULL);

    // ******* Copying Christoffel symbol and Ricci Tensor code here ********

    // Copying for compute_christoffel
    Tensor<3, data_t, 3> ULL;        //!< standard christoffel symbols
    Tensor<3, data_t, 3> LLL;        //!< 3 lower indices
    Tensor<1, data_t, 3> contracted; //!< contracted christoffel
    FOR_3D(i, j, k)
    {
        LLL[i][j][k] = 0.5 * (vars.d1_h[j][i][k] + vars.d1_h[k][i][j] -
                              vars.d1_h[j][k][i]);
    }
    FOR_3D(i, j, k)
    {
        ULL[i][j][k] = 0;
        FOR_3D(l) { ULL[i][j][k] += h_UU[i][l] * LLL[l][j][k]; }
    }
    FOR_3D(i)
    {
        contracted[i] = 0;
        FOR_3D(j, k) { contracted[i] += h_UU[j][k] * ULL[i][j][k]; }
    }

    // Copying for compute_ricci

    Tensor<2, data_t, 3> ricci_LL; // Ricci with two indices down

    data_t boxtildechi = 0;

    Tensor<2, data_t, 3> covdtilde2chi;
    FOR_3D(k, l)
    {
        covdtilde2chi[k][l] = vars.d2_chi[k][l];
        FOR_3D(m) { covdtilde2chi[k][l] -= ULL[m][k][l] * vars.d1_chi[m]; }
    }

    FOR_3D(k, l) { boxtildechi += h_UU[k][l] * covdtilde2chi[k][l]; }

    data_t dchi_dot_dchi = 0;
    {
        FOR_3D(m, n)
        {
            dchi_dot_dchi += h_UU[m][n] * vars.d1_chi[m] * vars.d1_chi[n];
        }
    }

    FOR_3D(i, j)
    {
        data_t ricci_tilde = 0;
        FOR_3D(k)
        {
            // Trick: For CCZ4, we can add Z terms to ricci by changing
            // Gamma to chrisvec This way of writing it allows the user to
            // pass Z/chi = {0};
            // Z_over_chi passed as 0 here
            ricci_tilde += 0.5 * (vars.h[k][i] * vars.d1_Gamma[k][j] +
                                  vars.h[k][j] * vars.d1_Gamma[k][i]);
            ricci_tilde +=
                0.5 * (vars.Gamma[k]) * (LLL[i][j][k] + LLL[j][i][k]);
            FOR_3D(l)
            {
                ricci_tilde -= 0.5 * h_UU[k][l] * vars.d2_h[i][j][k][l];
                FOR_3D(m)
                {
                    ricci_tilde += h_UU[l][m] * (ULL[k][l][i] * LLL[j][k][m] +
                                                 ULL[k][l][j] * LLL[i][k][m] +
                                                 ULL[k][i][m] * LLL[k][l][j]);
                }
            }
        }

        data_t ricci_chi =
            0.5 * ((3 - 2) * covdtilde2chi[i][j] + vars.h[i][j] * boxtildechi -
                   ((3 - 2) * vars.d1_chi[i] * vars.d1_chi[j] +
                    3 * vars.h[i][j] * dchi_dot_dchi) /
                       (2 * vars.chi));

        data_t z_terms = 0;
        /*FOR(k)
        {
            z_terms +=
                Z_over_chi[k] *
                (vars.h[i][k] * d1.chi[j] + vars.h[j][k] * d1.chi[i] -
                 vars.h[i][j] * d1.chi[k] + d1.h[i][j][k] * vars.chi);
        }*/

        ricci_LL[i][j] = (ricci_chi + vars.chi * ricci_tilde) /
                         vars.chi; //+ z_terms) / vars.chi;
    }

    // Copying for computing_phys_chris
    Tensor<3, data_t, 3> chris_phys;
    FOR_3D(i, j, k)
    {
        chris_phys[i][j][k] = ULL[i][j][k] - 0.5 / vars.chi *
                                                 (delta(i, k) * vars.d1_chi[j] +
                                                  delta(i, j) * vars.d1_chi[k]);
        FOR_3D(m)
        {
            chris_phys[i][j][k] +=
                0.5 / vars.chi * vars.h[j][k] * h_UU[i][m] * vars.d1_chi[m];
        }
    }

    // ******** Finished Copying *******

    // Extrinsic curvature and corresponding covariant and partial derivatives
    FOR_3D(i, j)
    {
        K_tensor[i][j] = vars.A[i][j] / vars.chi +
                         1. / 3. * (vars.h[i][j] * vars.K) / vars.chi;

        FOR_3D(k)
        {
            d1_K_tensor[i][j][k] =
                vars.d1_A[i][j][k] / vars.chi -
                vars.d1_chi[k] / vars.chi * K_tensor[i][j] +
                1. / 3. * vars.d1_h[i][j][k] * vars.K / vars.chi +
                1. / 3. * vars.h[i][j] * vars.d1_K[k] / vars.chi;
        }
    }
    // covariant derivative of K
    FOR_3D(i, j, k)
    {
        covariant_deriv_K_tensor[i][j][k] = d1_K_tensor[i][j][k];
    }

    FOR_3D(i, j, k, l)
    {
        covariant_deriv_K_tensor[i][j][k] +=
            -chris_phys[l][k][i] * K_tensor[l][j] -
            chris_phys[l][k][j] * K_tensor[i][l];
    }

    // Calculate electric and magnetic fields
    FOR_3D(i, j)
    {
        out.E[i][j] = 0;
        out.B[i][j] = 0;
    }

    FOR_3D(i, j, k, l)
    {
        out.B[i][j] +=
            epsilon3_LUU[i][k][l] * (covariant_deriv_K_tensor[l][j][k]);
    }

    FOR_3D(i, j) { out.E[i][j] += ricci_LL[i][j] + vars.K * K_tensor[i][j]; }

    FOR_3D(i, j, k, l)
    {
        out.E[i][j] += -K_tensor[i][k] * K_tensor[l][j] * h_UU[k][l] * vars.chi;
    }

    return out;
}

// Calculation of the Weyl4 scalar
template <class data_t>
NPScalar_t<data_t> Weyl4::compute_Weyl4(const EBFields_t<data_t> &ebfields,
                                        const CCZ4vars3D<data_t> &vars,
                                        const Coordinates<data_t> &coords) const
{
    NPScalar_t<data_t> out;

    // Calculate the tetrads
    const Tetrad_t<data_t> tetrad = compute_null_tetrad(vars, coords);

    // Projection of Electric and magnetic field components using tetrads
    out.Real = 0.0;
    out.Im = 0.0;
    FOR_3D(i, j)
    {
        out.Real += 0.5 * (ebfields.E[i][j] * (tetrad.w[i] * tetrad.w[j] -
                                               tetrad.v[i] * tetrad.v[j]) -
                           2.0 * ebfields.B[i][j] * tetrad.w[i] * tetrad.v[j]);
        out.Im += 0.5 * (ebfields.B[i][j] * (-tetrad.w[i] * tetrad.w[j] +
                                             tetrad.v[i] * tetrad.v[j]) -
                         2.0 * ebfields.E[i][j] * tetrad.w[i] * tetrad.v[j]);
    }

    return out;
}

// Calculation of the null tetrad
// Defintions from gr-qc/0104063
// "The Lazarus project: A pragmatic approach to binary black hole evolutions",
// Baker et al.
template <class data_t>
Tetrad_t<data_t>
Weyl4::compute_null_tetrad(const CCZ4vars3D<data_t> &vars,
                           const Coordinates<data_t> &coords) const
{
    Tetrad_t<data_t> out;

    // compute coords
    const data_t x = coords.x;
    const double y = coords.y;
    const double z = 0; // coords.z;

    // the inverse metric and alternating levi civita symbol
    data_t deth = vars.h[0][0] * vars.h[1][1] * vars.h[2][2] +
                  2 * vars.h[0][1] * vars.h[0][2] * vars.h[1][2] -
                  vars.h[0][0] * vars.h[1][2] * vars.h[1][2] -
                  vars.h[1][1] * vars.h[0][2] * vars.h[0][2] -
                  vars.h[2][2] * vars.h[0][1] * vars.h[0][1];
    data_t deth_inverse = 1. / deth;
    Tensor<2, data_t, 3> h_UU;
    h_UU[0][0] = (vars.h[1][1] * vars.h[2][2] - vars.h[1][2] * vars.h[1][2]) *
                 deth_inverse;
    h_UU[0][1] = (vars.h[0][2] * vars.h[1][2] - vars.h[0][1] * vars.h[2][2]) *
                 deth_inverse;
    h_UU[0][2] = (vars.h[0][1] * vars.h[1][2] - vars.h[0][2] * vars.h[1][1]) *
                 deth_inverse;
    h_UU[1][1] = (vars.h[0][0] * vars.h[2][2] - vars.h[0][2] * vars.h[0][2]) *
                 deth_inverse;
    h_UU[1][2] = (vars.h[0][1] * vars.h[0][2] - vars.h[0][0] * vars.h[1][2]) *
                 deth_inverse;
    h_UU[2][2] = (vars.h[0][0] * vars.h[1][1] - vars.h[0][1] * vars.h[0][1]) *
                 deth_inverse;
    h_UU[1][0] = h_UU[0][1];
    h_UU[2][0] = h_UU[0][2];
    h_UU[2][1] = h_UU[1][2];

    // Copy the epsilon function for general dimension here
    Tensor<3, double, 3> epsilon = {0.};
    epsilon[0][1][2] = 1.0;
    epsilon[1][2][0] = 1.0;
    epsilon[2][0][1] = 1.0;
    epsilon[0][2][1] = -1.0;
    epsilon[2][1][0] = -1.0;
    epsilon[1][0][2] = -1.0;

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

    FOR_3D(i, j, k, m)
    {
        out.w[i] += pow(chi, -0.5) * h_UU[i][j] * epsilon[j][k][m] * out.v[k] *
                    out.u[m];
    }

    // Gram Schmitt orthonormalisation
    // Choice of orthonormalisaion to avoid frame-dragging
    data_t omega_11 = 0.0;
    FOR_3D(i, j) { omega_11 += out.v[i] * out.v[j] * vars.h[i][j] / chi; }
    FOR_3D(i) { out.v[i] = out.v[i] / sqrt(omega_11); }

    data_t omega_12 = 0.0;
    FOR_3D(i, j) { omega_12 += out.v[i] * out.u[j] * vars.h[i][j] / chi; }
    FOR_3D(i) { out.u[i] += -omega_12 * out.v[i]; }

    data_t omega_22 = 0.0;
    FOR_3D(i, j) { omega_22 += out.u[i] * out.u[j] * vars.h[i][j] / chi; }
    FOR_3D(i) { out.u[i] = out.u[i] / sqrt(omega_22); }

    data_t omega_13 = 0.0;
    data_t omega_23 = 0.0;
    FOR_3D(i, j)
    {
        omega_13 += out.v[i] * out.w[j] * vars.h[i][j] / chi;
        omega_23 += out.u[i] * out.w[j] * vars.h[i][j] / chi;
    }
    FOR_3D(i) { out.w[i] += -(omega_13 * out.v[i] + omega_23 * out.u[i]); }

    data_t omega_33 = 0.0;
    FOR_3D(i, j) { omega_33 += out.w[i] * out.w[j] * vars.h[i][j] / chi; }
    FOR_3D(i) { out.w[i] = out.w[i] / sqrt(omega_33); }

    return out;
}

#endif /* WEYL4_HPP_ */
