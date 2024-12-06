/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CONSTRAINTSCARTOON_HPP_)
#error "This file should only be included through ConstraintsCartoon.hpp"
#endif

#ifndef CONSTRAINTSCARTOON_IMPL_HPP_
#define CONSTRAINTSCARTOON_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

template <class coupling_t>
inline Constraints<coupling_t>::Constraints(
    double dx, coupling_t a_coupling, double a_G_Newton,
    double cosmological_constant /*defaulted*/)
    : m_deriv(dx), m_coupling(a_coupling), m_G_Newton(a_G_Newton),
      m_cosmological_constant(cosmological_constant), m_dx(dx)
{
}

template <class coupling_t>
template <class data_t>
void Constraints<coupling_t>::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    Coordinates<data_t> coords{current_cell, m_dx};

    // Function below calls function to compute matter terms, too
    constraints_t<data_t> out = constraint_equations(vars, d1, d2, coords.y);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Ham, c_Ham);
#if CH_SPACEDIM == 3
    current_cell.store_vars(out.Mom, GRInterval<c_Mom1, c_Mom3>());
#elif CH_SPACEDIM == 2
    current_cell.store_vars(out.Mom, GRInterval<c_Mom1, c_Mom2>());
#else
#ifdef CH_SPACEDIM
#error current_cell.store_vars() has not got your dimension combination implemented.
#endif
#endif

    current_cell.store_vars(out.Mom_abs, c_Mom);
    current_cell.store_vars(out.rho, c_rho);
    current_cell.store_vars(out.rho_ADM, c_rho_ADM);
    current_cell.store_vars(out.Si[0], c_Sx);
    current_cell.store_vars(out.Si[1], c_Sy);
    // current_cell.store_vars(out.Sij[0][0], c_Sxx);
    // current_cell.store_vars(out.Sij[0][1], c_Sxy);
    // current_cell.store_vars(out.Sij[1][1], c_Syy);
    // current_cell.store_vars(out.Sww, c_Sww);
    // current_cell.store_vars(out.S, c_S);

    // current_cell.store_vars(out.Sij_TF[0][0], c_Sxx_TF);
    // current_cell.store_vars(out.Sij_TF[0][1], c_Sxy_TF);
    // current_cell.store_vars(out.Sij_TF[1][1], c_Syy_TF);
    // current_cell.store_vars(out.Sww_TF, c_Sww_TF);
}

template <class coupling_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
typename Constraints<coupling_t>:: template constraints_t<data_t>
Constraints<coupling_t>::constraint_equations(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const double &cartoon_coord) const
{
    using namespace TensorAlgebra;

    constraints_t<data_t> out;

    const int dI = CH_SPACEDIM - 1;
    const int nS =
        GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
    const double one_over_cartoon_coord = 1. / cartoon_coord;

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    auto h_UU_ww = 1. / vars.hww;
    auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    Tensor<1, data_t> chris_ww;
    FOR(i)
    {
        chris_ww[i] =
            one_over_cartoon_coord * (delta(i, dI) - h_UU[i][dI] * vars.hww);
        FOR(j) chris_ww[i] -= 0.5 * h_UU[i][j] * d1.hww[j];
    }

    Tensor<1, data_t> reg_07;
    FOR(i)
    {
        reg_07[i] =
            one_over_cartoon_coord * (vars.A[i][dI] - delta(i, dI) * vars.Aww);
    }

    auto ricci = CCZ4CartoonGeometry::compute_ricci(vars, d1, d2, h_UU, h_UU_ww,
                                                    chris, cartoon_coord);

    auto A_UU = TensorAlgebra::raise_all(vars.A, h_UU);
    data_t tr_A2 = TensorAlgebra::compute_trace(vars.A, A_UU) +
                   nS * vars.Aww * vars.Aww * h_UU_ww * h_UU_ww;

    out.Ham = ricci.scalar +
              (GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM - tr_A2;
    out.Ham -= 2 * m_cosmological_constant;

    Tensor<2, data_t> covd_A[CH_SPACEDIM];
    FOR(i, j, k)
    {
        covd_A[i][j][k] = d1.A[j][k][i];
        FOR(l)
        {
            covd_A[i][j][k] += -chris.ULL[l][i][j] * vars.A[l][k] -
                               chris.ULL[l][i][k] * vars.A[l][j];
        }
    }

    FOR(i)
    {
        out.Mom[i] =
            -(GR_SPACEDIM - 1.) * d1.K[i] / GR_SPACEDIM +
            nS * h_UU_ww * (reg_07[i] - 0.5 * h_UU_ww * vars.Aww * d1.hww[i]);
    }
    FOR(i, j)
    {
        out.Mom[i] -= nS * h_UU_ww * chris_ww[j] * vars.A[i][j];
        FOR(k)
        {
            out.Mom[i] += h_UU[j][k] * (covd_A[k][j][i] -
                                        GR_SPACEDIM * vars.A[i][j] * d1.chi[k] /
                                            (2 * chi_regularised));
        }
    }

    // Matter contributions - probably not needed
    // data_t coupling =0., f_coupling=0., fprime=0.;
    // m_coupling.compute_coupling(f_coupling, fprime, coupling, vars);

    emtensorCartoon_t<data_t> emtensor = compute_EMS_EM_tensor(vars, d1, h_UU,
                                               h_UU_ww, chris, nS, m_coupling);
    out.Ham += -16.0 * M_PI * m_G_Newton * emtensor.rho;
    FOR(i)
    {
        out.Mom[i] += -8.0 * M_PI * m_G_Newton * emtensor.Si[i];
    }

    data_t Mom_abs = 0.;

    FOR(i) { Mom_abs += out.Mom[i] * out.Mom[i]; }
    Mom_abs = sqrt(Mom_abs);
    out.Mom_abs = Mom_abs;

    out.rho = emtensor.rho;
    data_t det_gamma = TensorAlgebra::compute_determinant(vars.h);
    det_gamma *= vars.hww;
    out.rho_ADM = out.rho * det_gamma;
    FOR(i)
    {
        out.Si[i] = emtensor.Si[i];
    }

    return out;
}

#endif /* CONSTRAINTSCARTOON_IMPL_HPP_ */
