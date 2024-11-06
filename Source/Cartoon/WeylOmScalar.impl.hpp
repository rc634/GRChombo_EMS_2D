/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(WEYLOMSCALAR_HPP_)
#error "This file should only be included through WeylOmScalar.hpp"
#endif

#ifndef WEYLOMSCALAR_IMPL_HPP_
#define WEYLOMSCALAR_IMPL_HPP_

#include "simd.hpp"

template <class data_t>
void WeylOmScalar::compute(Cell<data_t> current_cell) const
{

    // copy data from chombo gridpoint into local variables
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

    // Get the coordinates
    const Coordinates<data_t> coords(current_cell, m_deriv.m_dx, m_center);

    // Compute the relevant components of the spacetime Riemann tensor
    Riemann_t<data_t> full_riemann =
        compute_riemann_tensor(vars, d1, d2, coords);

    // work out the Newman Penrose scalar matrix
    NPScalarMatrix_t<data_t> out =
        compute_Weyl_Om(full_riemann, vars, d1, d2, coords);

    // Debugging
    /*if(abs(coords.x-0.53)<m_deriv.m_dx && abs(coords.y-0.41)<m_deriv.m_dx &&
    abs(coords.z-0.75)<m_deriv.m_dx)
    {
           std::cout.precision(15);
           std::cout << " x = "<< coords.x << " y = "<< coords.y << " z = "<<
    coords.z << " Om22 = " << out.Om22 << " Om23 = " << out.Om23 << " Om33 = "
    << out.Om33 << " Omww = " << out.Omww << std::endl;
    }*/

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Om22, c_Weyl4_Re);
    current_cell.store_vars(0., c_Weyl4_Im);
    //    current_cell.store_vars(out.Om22[1], c_Weyl4_Im);
    //    current_cell.store_vars(out.Om33, c_Weyl_Om33);
    //    current_cell.store_vars(out.Omww, c_Weyl_Omww);
}

// Calculation of the Weyl scalar matrix
template <class data_t>
NPScalarMatrix_t<data_t> WeylOmScalar::compute_Weyl_Om(
    const Riemann_t<data_t> &full_riemann, const Vars<data_t> &vars,
    const Vars<Tensor<1, data_t>> &d1, const Diff2Vars<Tensor<2, data_t>> &d2,
    const Coordinates<data_t> &coords) const
{
    NPScalarMatrix_t<data_t> out;

    // Calculate null tetrads
    const Tetrad_t<data_t> tetrad = compute_null_tetrad(vars, coords);

    Tensor<2, data_t> mvec = 0.;
    FOR(i)
    {
        mvec[0][i] = tetrad.m1[i];
        mvec[1][i] = tetrad.m2[i];
        //        mvec[2][i] = tetrad.m3[i];
        //	cout<<mvec[0][i]<<endl;
        //	cout<<mvec[1][i]<<endl;
    }

    Tensor<2, data_t> omlocij = 0.;
    for (int i = 1; i < CH_SPACEDIM; i++)
    {
        for (int j = 1; j < CH_SPACEDIM; j++)
        {
            FOR(k, l)
            {
                omlocij[i][j] +=
                    full_riemann.riemann_0d0d[k][l] * mvec[i][k] * mvec[j][l];
		//std::cout << "0d0d   " << omlocij[i][j] << "\n" << std::endl;
                FOR(m)
                {
                    omlocij[i][j] += full_riemann.riemann_0ddd[l][k][m] *
                                         mvec[j][l] * mvec[i][k] * mvec[0][m] +
                                     full_riemann.riemann_0ddd[k][l][m] *
                                         mvec[i][k] * mvec[j][l] * mvec[0][m];
 		    //std::cout << "0d0d   " << omlocij[i][j] << "\n" << std::endl;
                    FOR(n)
                    {
                        omlocij[i][j] += full_riemann.riemann_dddd[m][k][n][l] *
                                         mvec[0][m] * mvec[i][k] * mvec[0][n] *
                     			 mvec[j][l];
			//std::cout << "dddd   " << omlocij[i][j] << "\n" << std::endl;
                    }
                }
            }
            omlocij[i][j] *= 0.25;
        }
    }

    // Compute the non-zero components of the Weyl Scalar Matrix: Eqs. (4.55),
    // (4.57)
    out.Om22 = omlocij[1][1];
    //    out.Om23 = omlocij[1][2];
    //    out.Om33 = omlocij[2][2];

    //    cout<<out.Om22<<endl;

    data_t gww = vars.hww / vars.chi;
    out.Omww = full_riemann.riemann_0w0w;
    FOR(k)
    {
        out.Omww -= 2. * full_riemann.riemann_0wdw[k] * tetrad.m1[k];
        FOR(l)
        {
            out.Omww +=
                full_riemann.riemann_dwdw[k][l] * tetrad.m1[k] * tetrad.m1[l];
        }
    }
    out.Omww = 0.25 * out.Omww / gww;

    out.Om22 -= out.Omww;

    // Debugging
    /*if(abs(coords.x-0.53)<m_deriv.m_dx && abs(coords.y-0.41)<m_deriv.m_dx &&
    abs(coords.z-0.75)<m_deriv.m_dx)
    {
        std::cout.precision(15);
        std::cout << " x = "<< coords.x << " y = "<< coords.y << " z = "<<
    coords.z << " Om22 = " << out.Om22 << " Om23 = " << out.Om23 << " Om33 = "
    << out.Om33 << " Omww = " << out.Omww << std::endl;
    }*/

    return out;
}

// Calculation of the null tetrad
// Defintions from arXiv:1609.01292
// "Extraction of gravitational-wave energy in higher dimensional numerical
// relativity using the Weyl tensor ", Cook and Sperhake
template <class data_t>
Tetrad_t<data_t>
WeylOmScalar::compute_null_tetrad(const Vars<data_t> &vars,
                                  const Coordinates<data_t> &coords) const
{
    Tetrad_t<data_t> out;

    // compute coords
    const data_t x = coords.x;
    const double y = 0.;
    const double z = coords.y;

    // calculate the tetrad
    // tetrad adapted to the standard highe dimensional sphere
    /*out.m1[0] = x;
    out.m1[1] = y;
    out.m1[2] = z;

    out.m2[0] = -y*y - z*z;
    out.m2[1] = x*y;
    out.m2[2] = x*z;

    out.m3[0] = 0.;
    out.m3[1] = -z;
    out.m3[2] = y;*/

    // tetrad as in the 4d code
    out.m1[0] = x;
    //    out.m1[1] = y;
    out.m1[1] = z;

    //    out.m2[0] = x*z;
    //    out.m2[1] = y*z;
    //    out.m2[1] = -x*x -y*y;

    out.m2[0] = -z;
    out.m2[1] = x;

    out.m3[0] = -y;
    //    out.m3[1] = x;
    out.m3[1] = 0.;

    // floor on chi
    const data_t chi = simd_max(vars.chi, 1e-4);
    // spacetime metric
    Tensor<2, data_t> gd = 0.;
    FOR(i, j) { gd[i][j] = vars.h[i][j] / chi; }
    data_t gww = vars.hww / chi;

    // Gram Schmitt orthonormalisation
    // Choice of orthonormalisaion to avoid frame-dragging
    data_t omega_11 = 0.;
    FOR(i, j) { omega_11 += out.m1[i] * out.m1[j] * gd[i][j]; }
    FOR(i) { out.m1[i] = out.m1[i] / sqrt(omega_11); }

    data_t omega_12 = 0.;
    FOR(i, j) { omega_12 += out.m2[i] * out.m1[j] * gd[i][j]; }
    FOR(i) { out.m2[i] += -omega_12 * out.m1[i]; }

    data_t omega_22 = 0.;
    FOR(i, j) { omega_22 += out.m2[i] * out.m2[j] * gd[i][j]; }
    FOR(i) { out.m2[i] = out.m2[i] / sqrt(omega_22); }

    data_t omega_13 = 0.;
    data_t omega_23 = 0.;
    FOR(i, j)
    {
        omega_13 += out.m3[i] * out.m1[j] * gd[i][j];
        omega_23 += out.m3[i] * out.m2[j] * gd[i][j];
    }
    FOR(i) { out.m3[i] += -(omega_13 * out.m1[i] + omega_23 * out.m2[i]); }

    data_t omega_33 = 0.;
    FOR(i, j) { omega_33 += out.m3[i] * out.m3[j] * gd[i][j]; }
    FOR(i) { out.m3[i] = out.m3[i] / sqrt(omega_33); }

    // The extra-dimensional component is already orthogonal
    out.mw = sqrt(1. / gww);

    return out;
}

// Calculation of the spatial components of the Riemann tensor
// with respect to the physical metric as in arXiv:1609.01292
template <class data_t>
Riemann_t<data_t>
WeylOmScalar::compute_riemann_tensor(const Vars<data_t> &vars,
                                     const Vars<Tensor<1, data_t>> &d1,
                                     const Diff2Vars<Tensor<2, data_t>> &d2,
                                     const Coordinates<data_t> &coords) const
{
    Riemann_t<data_t> out;

    // const int nS =
    //     GR_SPACEDIM -
    //     CH_SPACEDIM /*2 debugging*/; //!< Dimensions of the transverse sphere
    const int dims =
        GR_SPACEDIM + 1 /*6 debugging*/; //!< Total Spacetime dimensions

    // compute position on grid relative to center
    // const data_t x = coords.x;
    // const double y = 0.;        // coords.y;
    const double z = coords.y; // coords.z;

    const double z2 = z * z;
    const double one_over_z = 1. / z;
    const double one_over_z2 = 1. / (z * z);

    // the inverse conformal metric and the corresponding Christoffel symbols
    using namespace TensorAlgebra;
    // const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    // const auto chris_t = compute_christoffel(d1.h, h_UU);

    // physical metric and its derivatives
    data_t gww, guww, Kww;
    Tensor<1, data_t> d1_gww, d1_Kww;
    Tensor<2, data_t> gd, Kdd, d2_gww;
    Tensor<2, Tensor<1, data_t>> d1_gd, d1_Kdd, DKdd;
    Tensor<4, data_t> d2_gd;

    const auto one_over_chi = 1. / simd_max(vars.chi, 1e-8);
    const auto one_over_chi2 = 1. / pow(vars.chi, 2);
    const auto one_over_chi3 = 1. / pow(vars.chi, 3);
    const double one_over_gr_spacedim =
        1. / ((double)GR_SPACEDIM) /*1./(dims - 1.) Debugging*/;

    gww = one_over_chi * vars.hww;
    guww = 1. / gww;
    Kww = one_over_chi * (vars.Aww + one_over_gr_spacedim * vars.hww * vars.K);

    FOR(i, j)
    {
        gd[i][j] = one_over_chi * vars.h[i][j];
        Kdd[i][j] = one_over_chi * (vars.A[i][j] + one_over_gr_spacedim *
                                                       vars.h[i][j] * vars.K);
    }

    FOR(i)
    {
        d1_gww[i] =
            one_over_chi * d1.hww[i] - one_over_chi2 * vars.hww * d1.chi[i];
        d1_Kww[i] =
            one_over_chi *
            (d1.Aww[i] +
             one_over_gr_spacedim * (vars.hww * d1.K[i] + vars.K * d1.hww[i]) -
             Kww * d1.chi[i]);
    }

    FOR(i, j)
    {
        d2_gww[i][j] =
            one_over_chi * d2.hww[i][j] -
            one_over_chi2 * (d1.hww[i] * d1.chi[j] + d1.hww[j] * d1.chi[i] +
                             vars.hww * d2.chi[i][j]) +
            2.0 * one_over_chi3 * vars.hww * d1.chi[i] * d1.chi[j];
        FOR(k)
        {
            d1_gd[i][j][k] = one_over_chi * d1.h[i][j][k] -
                             one_over_chi2 * vars.h[i][j] * d1.chi[k];
            d1_Kdd[i][j][k] = one_over_chi *
                              (d1.A[i][j][k] +
                               one_over_gr_spacedim * (vars.h[i][j] * d1.K[k] +
                                                       vars.K * d1.h[i][j][k]) -
                               Kdd[i][j] * d1.chi[k]);
            FOR(l)
            {
                d2_gd[i][j][k][l] =
                    one_over_chi * d2.h[i][j][k][l] -
                    one_over_chi2 *
                        (d1.h[i][j][k] * d1.chi[l] + d1.h[i][j][l] * d1.chi[k] +
                         vars.h[i][j] * d2.chi[k][l]) +
                    2. * one_over_chi3 * vars.h[i][j] * d1.chi[k] * d1.chi[l];
            }
        }
    }

    const auto gu = TensorAlgebra::compute_inverse_sym(gd);
    const auto chris = compute_christoffel(d1_gd, gu);

    // Regularisation terms
    // Since we have one less spatial dimension,
    // we lower the rank for all regularisation term,
    // and changed all index 2 -> 1

    Tensor<1, data_t> reg07 = 0.;
    reg07[0] = Kdd[0][1] / z;
    reg07[1] = (Kdd[1][1] - Kww) / z;

    Tensor<1, data_t> reg02 = 0.;
    reg02[0] = -gu[0][1] * gww / z;
    reg02[1] = (1.0 - gu[1][1] * gww) / z;

    Tensor<2, data_t> reg13 = 0.;
    reg13[0][0] = d1_gd[0][1][0] / z;
    reg13[0][1] = 0.5 * ((d1_gd[1][1][0] + d1_gd[0][1][1] - d1_gww[0]) / z -
                         gd[0][1] / z2);
    reg13[1][1] = (d1_gd[1][1][1] - d1_gww[1]) / z - (gd[1][1] - gww) / z2;
    reg13[1][0] = reg13[0][1];

    Tensor<2, data_t> reg08 = 0.;
    FOR(i, j) { reg08[i][j] = -0.5 * d1_gd[i][j][1] / z; }
    FOR(k)
    {
        reg08[1][k] += 0.5 * gd[1][k] / z2;
        reg08[k][1] += 0.5 * gd[k][1] / z2;
    }
    reg08[1][1] -= gww / z2;

    data_t reg14 = 0.;
    FOR(i) { reg14 -= gww * gu[1][i] * d1_gww[i] * one_over_z; }
    data_t reg15 = one_over_z2 * (gww - gu[1][1] * gww * gww);

    // Christoffel symbols
    // Eq. (4.12)
    Tensor<1, data_t> chris_Uww = reg02;
    FOR(i)
    {
        FOR(j) { chris_Uww[i] -= 0.5 * gu[i][j] * d1_gww[j]; }
    }

    // covariant derivative of \nabla_i K_{jk}, with indices down
    // we follow the usual convention and the derivative index is the last one:
    // K_{jk;i}
    FOR(i, j, k)
    {
        DKdd[j][k][i] = d1_Kdd[j][k][i];
        FOR(l)
        {
            DKdd[j][k][i] += -chris.ULL[l][i][j] * Kdd[k][l] -
                             chris.ULL[l][i][k] * Kdd[j][l];
        }
    }

    // ---------------- Riemann3 ---------------------
    Tensor<4, data_t> R3ijkl = 0.;
    FOR(i, j, k, l)
    {
        R3ijkl[i][j][k][l] = 0.5 * (d2_gd[j][k][i][l] + d2_gd[i][l][j][k] -
                                    d2_gd[j][l][i][k] - d2_gd[i][k][j][l]);
        FOR(m, n)
        {
            R3ijkl[i][j][k][l] +=
                gd[m][n] * (chris.ULL[n][i][l] * chris.ULL[m][j][k] -
                            chris.ULL[n][i][k] * chris.ULL[m][j][l]);
        }
    }

    Tensor<2, data_t> R3wjwl = 0.;
    FOR(i, k)
    {
        R3wjwl[i][k] = reg13[i][k] + reg08[i][k] - 0.5 * d2_gww[i][k] +
                       0.25 * guww * d1_gww[i] * d1_gww[k];
        FOR(m, n)
        {
            R3wjwl[i][k] -= gd[m][n] * chris.ULL[n][i][k] * chris_Uww[m];
        }
    }

    data_t R3wuwu = reg14 + reg15;
    FOR(m, n) { R3wuwu -= 0.25 * gu[m][n] * d1_gww[m] * d1_gww[n]; }
    // -----------------------------------------------

    // ---------------- Ricci3 -----------------------
    // We are working on SO(D-2), so we D-d-1 = d-3
    Tensor<2, data_t> R3ij = 0.;
    FOR(i, j)
    {
        R3ij[i][j] = (dims - 3) * guww * R3wjwl[i][j];
        FOR(m, n) { R3ij[i][j] += gu[m][n] * R3ijkl[m][i][n][j]; }
    }

    data_t R3ww = (dims - 4) * guww * R3wuwu;
    FOR(m, n) { R3ww += gu[m][n] * R3wjwl[m][n]; }

    // Debugging
    /*if(abs(coords.x-0.53)<m_deriv.m_dx && abs(coords.y-0.41)<m_deriv.m_dx &&
    abs(coords.z-0.75)<m_deriv.m_dx)
    {
        std::cout.precision(17);
        std::cout << " dims = " << dims << std::endl;
        //std::cout << " x = "<< x << " y = "<< y << " z = "<< z << " R3wjwl = "
    << R3wjwl[0][0] << " " << R3wjwl[0][1] << " " << R3wjwl[0][2] << " " <<
    R3wjwl[1][1] << " " << R3wjwl[1][2] << " " << R3wjwl[2][2] << std::endl;
        //std::cout << " x = "<< x << " y = "<< y << " z = "<< z << " R3ijkl = "
    << R3ijkl[0][1][0][1] << " " << R3ijkl[0][1][0][2] << " " <<
    R3ijkl[0][1][1][2] << " " << R3ijkl[0][2][0][2] << " " << R3ijkl[0][2][1][2]
    << " " << R3ijkl[1][2][1][2] << std::endl;
        //std::cout << " x = "<< x << " y = "<< y << " z = "<< z << " R3ij = "
    << R3ij[0][0] << " " << R3ij[0][1] << " " << R3ij[0][2] << " " << R3ij[1][1]
    << " " << R3ij[1][2] << " " << R3ij[2][2] << std::endl; std::cout << " x =
    "<< x << " y = "<< y << " z = "<< z << " R3wuwu = " << R3wuwu <<  " R3ww = "
    << R3ww << std::endl;
    }*/
    // -----------------------------------------------

    // ---------------- Riemann4 ---------------------
    // data_t R4w0w0 = R3ww + (vars.K - guww * Kww) * Kww;
    out.riemann_0w0w = R3ww + (vars.K - guww * Kww) * Kww;

    FOR(i, k)
    {
        out.riemann_0d0d[i][k] = R3ij[i][k] + vars.K * Kdd[i][k];
        FOR(m, n)
        {
            out.riemann_0d0d[i][k] -= gu[m][n] * Kdd[i][m] * Kdd[k][n];
        }
    }

    Tensor<1, data_t> R4w0wl = 0.;
    FOR(i)
    {
        R4w0wl[i] = d1_Kww[i] - 0.5 * guww * Kww * d1_gww[i] - reg07[i];
        FOR(m) { R4w0wl[i] += chris_Uww[m] * Kdd[m][i]; }
        out.riemann_0wdw[i] = R4w0wl[i];
    }

    // Tensor<2, data_t> R4wjwl = 0.;
    FOR(i, j)
    {
        // R4wjwl[i][j] = R3wjwl[i][j] + Kdd[i][j] * Kww;
        out.riemann_dwdw[i][j] = R3wjwl[i][j] + Kdd[i][j] * Kww;
    }

    // Tensor<3, data_t> R4i0kl = 0.;
    FOR(i, k, l)
    {
        // R4i0kl[i][k][l] = DKdd[i][k][l] - DKdd[i][l][k];
        // out.riemann_0ddd[i][k][l] = - R4i0kl[i][k][l];
        out.riemann_0ddd[i][k][l] = DKdd[i][l][k] - DKdd[i][k][l];
    }

    // Tensor<4, data_t> R4ijkl = 0.;
    FOR(i, j, k, l)
    {
        // R4ijkl[i][j][k][l] = R3ijkl[i][j][k][l] + Kdd[i][k]*Kdd[j][l] -
        // Kdd[i][l]*Kdd[j][k];
        out.riemann_dddd[i][j][k][l] =
            R3ijkl[i][j][k][l] + Kdd[i][k] * Kdd[j][l] - Kdd[i][l] * Kdd[j][k];
    }

    // Debugging
    /*if(abs(coords.x-0.53)<m_deriv.m_dx && abs(coords.y-0.41)<m_deriv.m_dx &&
    abs(coords.z-0.75)<m_deriv.m_dx)
    {
        std::cout.precision(17);
        //std::cout << " x = "<< x << " y = "<< y << " z = "<< z << " R4i0k0 = "
    << out.riemann_0d0d[0][0] << " " << out.riemann_0d0d[0][1] << " " <<
    out.riemann_0d0d[0][2] << " " << out.riemann_0d0d[1][1] << " " <<
    out.riemann_0d0d[1][2] << " " << out.riemann_0d0d[2][2] << std::endl;
        std::cout << " x = "<< x << " y = "<< y << " z = "<< z << " R4iwkw = "
    << out.riemann_dwdw[0][0] << " " << out.riemann_dwdw[0][1] << " " <<
    out.riemann_dwdw[0][2] << " " << out.riemann_dwdw[1][1] << " " <<
    out.riemann_dwdw[1][2] << " " << out.riemann_dwdw[2][2] << std::endl;
        //std::cout << " x = "<< x << " y = "<< y << " z = "<< z << " R4w0wl = "
    << out.riemann_0wdw[0] << " " << out.riemann_0wdw[1] << " " <<
    out.riemann_0wdw[2] << std::endl;
        //std::cout << " x = "<< x << " y = "<< y << " z = "<< z << " R40ikl "
    << out.riemann_0ddd[2][0][1] << " " << out.riemann_0ddd[2][0][2] << " " <<
    out.riemann_0ddd[2][1][2] << std::endl;
        //std::cout << " x = "<< x << " y = "<< y << " z = "<< z << " R4ijkl "
    << out.riemann_dddd[1][2][0][1] << " " << out.riemann_dddd[1][2][0][2] << "
    " << out.riemann_dddd[1][2][1][2] << std::endl;
    }*/
    // NOTE: computation of Riemann4 agrees with Uli's code up to 10^-12

    // ---------------------------------------------------------

    return out;
}

#endif /* WEYLOMSCALAR_IMPL_HPP_ */
