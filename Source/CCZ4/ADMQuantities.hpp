/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMQUANTITIES_HPP_
#define ADMQUANTITIES_HPP_

#include "ADMConformalVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

#if CH_SPACEDIM == 2
#include "CCZ4CartoonVars.hpp"
#endif

//! Calculates the ADM mass and momentum (for 3D)
class ADMQuantities
{
    // Use the variable definition in ADMConformalVars - only require the key
    // vars
#if CH_SPACEDIM == 3
    template <class data_t> using Vars = ADMConformalVars::VarsNoGauge<data_t>;
    template <class data_t>
    using Diff1Vars = ADMConformalVars::Diff2VarsNoGauge<data_t>;
#elif CH_SPACEDIM == 2
    template <class data_t> using Vars = CCZ4CartoonVars::VarsNoGauge<data_t>;
    template <class data_t>
    using Diff1Vars = CCZ4CartoonVars::Diff2VarsNoGauge<data_t>;
#endif

  public:
    enum DIR
    {
        X,
        Y,
        Z
    };

#if CH_SPACEDIM == 3
    ADMQuantities(const std::array<double, CH_SPACEDIM> &a_center, double a_dx,
                  int a_c_Madm = -1, int a_c_Jadm = -1,
                  const Interval &a_c_Padm = Interval(),
                  double a_G_Newton = 1.0)
        : m_deriv(a_dx), m_center(a_center), m_G_Newton(a_G_Newton),
          m_c_Madm(a_c_Madm), m_c_Padm(a_c_Padm), m_c_Jadm(a_c_Jadm), m_dir(Z)
    {
        CH_assert(m_c_Padm.size() == 0 || m_c_Padm.size() == 1 ||
                  m_c_Padm.size() == GR_SPACEDIM);
    }
#elif CH_SPACEDIM == 2
    ADMQuantities(const std::array<double, CH_SPACEDIM> &a_center, double a_dx,
                  int a_c_Madm = -1, const Interval &a_c_Padm = Interval(),
                  double a_G_Newton = 1.0)
        : m_deriv(a_dx), m_center(a_center), m_G_Newton(a_G_Newton),
          m_c_Madm(a_c_Madm), m_c_Padm(a_c_Padm)
    {
        // in 2D, only compute norm of P, as Py=Pz=0
        CH_assert(m_c_Padm.size() == 0 || m_c_Padm.size() == 1);
    }
#endif

#if CH_SPACEDIM == 3
    // in case user wants to change direction of spin calculation to something
    // other than Z
    void set_spin_dir(DIR spin_direction) { m_dir = spin_direction; }
#endif

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and calc 1st
        // derivs
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Diff1Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);

        // Surface element for integration
        Coordinates<data_t> coords(current_cell, m_deriv.m_dx, m_center);
#if CH_SPACEDIM == 3
        Tensor<1, data_t> x = {coords.x, coords.y, coords.z};
#elif CH_SPACEDIM == 2
        Tensor<1, data_t> x = {coords.x, coords.y};
#endif
        Tensor<1, data_t> dS_U = x;

        data_t dS_norm = 0.;
        FOR(i, j) { dS_norm += vars.h[i][j] / vars.chi * dS_U[i] * dS_U[j]; }
        dS_norm = sqrt(dS_norm);
        FOR(i) { dS_U[i] /= dS_norm; }

        Tensor<1, data_t> dS_L;
        FOR(i)
        {
            // dS_L[i] = dS_U[i];
            dS_L[i] = 0.;
            FOR(j) { dS_L[i] += vars.h[i][j] / vars.chi * dS_U[j]; }
        }

        if (m_c_Madm >= 0)
        {
            data_t Madm = 0.0;
            data_t chi_3_over_2 = vars.chi * sqrt(vars.chi);
            FOR(i, j, k, l)
            {
                Madm += dS_L[i] / (16. * M_PI * m_G_Newton) / chi_3_over_2 *
                        h_UU[j][k] * h_UU[i][l] *
                        (vars.chi * (d1.h[l][k][j] - d1.h[j][k][l]) -
                         (vars.h[l][k] * d1.chi[j] - vars.h[j][k] * d1.chi[l]));
            }
#if CH_SPACEDIM == 2
            data_t h_UU_ww = 1. / vars.hww;
            const int dI = CH_SPACEDIM - 1;
            const double one_over_cartoon_coord = 1. / coords.y;
            FOR(l)
            {
                data_t dw_h_lw =
                    one_over_cartoon_coord *
                    (vars.h[l][dI] - TensorAlgebra::delta(l, dI) * vars.hww);
                FOR(i)
                {
                    Madm += dS_L[i] / (16. * M_PI * m_G_Newton) / chi_3_over_2 *
                            h_UU_ww * h_UU[i][l] *
                            (vars.chi * (dw_h_lw - d1.hww[l]) -
                             (-vars.hww * d1.chi[l]));
                }
            }
#endif

            // assign values of ADM Mass in output box
            current_cell.store_vars(Madm, m_c_Madm);
        }

        if (m_c_Padm.size() > 0)
        {
#if CH_SPACEDIM == 3
            Tensor<1, data_t> Padm = {0.};

            FOR(i, j)
            {
                Padm[i] += -dS_L[j] / (8. * M_PI * m_G_Newton) * vars.K *
                           TensorAlgebra::delta(i, j);

                FOR(k)
                {
                    Padm[i] += dS_L[j] / (8. * M_PI * m_G_Newton) * h_UU[j][k] *
                               (vars.A[k][i] + vars.K * vars.h[k][i] / 3.);
                }
            }

            // assign values of ADM Linear Momentum in output box
            if (m_c_Padm.size() == GR_SPACEDIM)
            {
                FOR(i)
                {
                    int ivar = m_c_Padm.begin() + i;
                    current_cell.store_vars(Padm[i], ivar);
                }
            }
            else if (m_c_Padm.size() == 1)
            {
                data_t Padm_sq = 0.0;
                FOR(i) { Padm_sq += Padm[i] * Padm[i]; }
                data_t Padm_abs = sqrt(Padm_sq);
                current_cell.store_vars(Padm_abs, m_c_Padm.begin());
            }
#elif CH_SPACEDIM == 2
            data_t Px_adm = 0.;

            int i = 0; // dir = x
            FOR(j)
            {
                Px_adm += -dS_L[j] / (8. * M_PI * m_G_Newton) * vars.K *
                          TensorAlgebra::delta(i, j);

                FOR(k)
                {
                    Px_adm += dS_L[j] / (8. * M_PI * m_G_Newton) * h_UU[j][k] *
                              (vars.A[k][i] + vars.K * vars.h[k][i] / 3.);
                }
            }

            // assign values of ADM Linear Momentum in output box
            current_cell.store_vars(Px_adm, m_c_Padm.begin());
#endif
        }

#if CH_SPACEDIM == 3
        if (m_c_Jadm >= 0)
        {
            // spin about m_dir axis (x, y or z)
            data_t Jadm = 0.0;

            // note this is the levi civita symbol,
            // not tensor (eps_tensor = eps_symbol * chi^-1.5)
            const Tensor<3, double> epsilon = TensorAlgebra::epsilon();

            FOR(i, j, k)
            {
                Jadm += -dS_L[i] / (8. * M_PI * m_G_Newton) *
                        epsilon[m_dir][j][k] * x[j] * vars.K *
                        TensorAlgebra::delta(i, k);

                FOR(l, m)
                {
                    Jadm += dS_L[i] / (8. * M_PI * m_G_Newton) *
                            epsilon[m_dir][j][k] * x[j] * h_UU[i][l] *
                            h_UU[k][m] * vars.chi *
                            (vars.A[l][m] + vars.K * vars.h[l][m] / 3.);
                }
            }

            // assign values of ADM Momentum in output box
            current_cell.store_vars(Jadm, m_c_Jadm);
        }
#endif
    }

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const std::array<double, CH_SPACEDIM> &m_center;
    const double m_G_Newton; //!< Newton's constant
    const int m_c_Madm;
    const Interval m_c_Padm;

#if CH_SPACEDIM == 3
    const int m_c_Jadm;
    DIR m_dir;
#endif
};

#endif /* ADMQUANTITIES_HPP_ */
