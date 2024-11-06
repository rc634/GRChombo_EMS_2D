/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This file calculates CCZ4 geometric quantities including cartoon terms (or a
// similar 3+1 split).
#ifndef CCZ4CARTOONGEOMETRY_HPP_
#define CCZ4CARTOONGEOMETRY_HPP_

#include "CCZ4Geometry.hpp"
#include "DimensionDefinitions.hpp"

//! A structure for the decomposed elements of the Energy Momentum Tensor in
//! 3+1D
template <class data_t> struct emtensorCartoon_t
{
    Tensor<2, data_t> Sij; //!< S_ij = T_ij
    Tensor<1, data_t> Si;  //!< S_i = T_ia_n^a
    data_t S;              //!< S = S^i_i
    data_t rho;            //!< rho = T_ab n^a n^b
    data_t Sww;            //!< S_ww = T^ww
};

template <class data_t> struct sf_potential
{
    data_t coupling;
    data_t dcoupling_dphi;
};

template <class data_t> struct ricciCartoon_t
{
    Tensor<2, data_t> LL; // Ricci with two indices down
    data_t scalar;        // Ricci scalar
    data_t LLww;          //!< ricci_ww with indices down
};

class CCZ4CartoonGeometry : CCZ4Geometry
{
  public:
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    static ricciCartoon_t<data_t> compute_ricci_Z(
        const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
        const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU, const data_t &h_UU_ww,
        const chris_t<data_t> &chris, const Tensor<1, data_t> &Z_over_chi,
        // const double &m_dx, //Debugging
        // const std::array<double, CH_SPACEDIM> &m_center, //Debugging
        // const data_t &x,//Debugging
        // const double &y,//Debugging
        const double &cartoon_coord)
    {
        using namespace TensorAlgebra;
        // Get standard 3+1 Ricci
        // Note that Z_over_chi includes the extra-dimensional terms of the
        // contracted Christoffels
        // auto ricci = CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris,
        // Z_over_chi);

        const int d = CH_SPACEDIM;
        const int dI = CH_SPACEDIM - 1;
        const int nS =
            GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
                                       //!< // Add cartoon terms to Ricci
        const double one_over_cartoon_coord = 1. / cartoon_coord;
        const double one_over_cartoon_coord2 =
            1. / (cartoon_coord * cartoon_coord);

        // NOTE: in the new code code we compute Ricci whilst in the old one we
        // computed Ricci * chi
        ricciCartoon_t<data_t> out;

        // regularisation terms
        Tensor<2, data_t> reg_08;
        FOR(i, j)
        {
            reg_08[i][j] = -0.5 * one_over_cartoon_coord * d1.h[i][j][dI] +
                           one_over_cartoon_coord2 *
                               (0.5 * (delta(i, dI) * vars.h[j][dI] +
                                       delta(j, dI) * vars.h[i][dI]) -
                                delta(i, dI) * delta(j, dI) * vars.hww);
        }

        Tensor<2, data_t> reg_09;
        FOR(i, j)
        {
            reg_09[i][j] =
                0.5 * one_over_cartoon_coord *
                ((h_UU_ww * vars.h[i][dI] - delta(i, dI)) * d1.hww[j] +
                 (h_UU_ww * vars.h[j][dI] - delta(j, dI)) * d1.hww[i]);
        }

        // Debugging - starts here. Alternative computation of Ricci as in the
        // old code
        Tensor<1, data_t> chris_ww;
        FOR(i)
        {
            chris_ww[i] = one_over_cartoon_coord *
                          (delta(i, dI) - h_UU[i][dI] * vars.hww);
            FOR(j) chris_ww[i] -= 0.5 * h_UU[i][j] * d1.hww[j];
        }
        Tensor<1, data_t>
            chris_contracted; //!< includes the higher D contributions
        FOR(i)
        {
            chris_contracted[i] =
                chris.contracted[i] + nS * h_UU_ww * chris_ww[i];
        }

        Tensor<2, data_t> c_ri_hh;
        FOR(i, j)
        {
            c_ri_hh[i][j] = nS * h_UU_ww *
                            (reg_08[i][j] + reg_09[i][j] -
                             0.25 * h_UU_ww * d1.hww[i] * d1.hww[j]);
            FOR(k)
            {
                c_ri_hh[i][j] +=
                    0.5 * (vars.h[i][k] * d1.Gamma[k][j] +
                           vars.h[j][k] * d1.Gamma[k][i] +
                           chris_contracted[k] *
                               (chris.LLL[i][j][k] + chris.LLL[j][i][k]));
                FOR(l)
                {
                    c_ri_hh[i][j] -= 0.5 * h_UU[k][l] * d2.h[i][j][k][l];
                    FOR(m)
                    {
                        c_ri_hh[i][j] +=
                            h_UU[l][m] *
                            (chris.ULL[k][l][i] * chris.LLL[j][k][m] +
                             chris.ULL[k][l][j] * chris.LLL[i][k][m] +
                             chris.ULL[k][i][m] * chris.LLL[k][l][j]);
                    }
                }
            }
        }
        // Debugging ends here

        data_t hUU_dchi_cartoon = 0;
        FOR(k) { hUU_dchi_cartoon += h_UU[dI][k] * d1.chi[k]; }
        hUU_dchi_cartoon = one_over_cartoon_coord * hUU_dchi_cartoon;

        // data_t dhww_dot_dchi = 0;
        // {
        //     FOR(m, n) { dhww_dot_dchi += h_UU[m][n] * d1.hww[m] * d1.chi[n];
        //     }
        // }

        // Add extra-dimensional terms to Ricci
        /*FOR(i, j)
        {
            data_t ricci_tilde_extra = 0;
            data_t ricci_chi_extra = 0;

            ricci_tilde_extra = h_UU_ww *
                                ( reg_08[i][j] + reg_09[i][j]
                                 -0.25 * h_UU_ww * d1.hww[i] * d1.hww[j]); //
        new cartoon terms

            ricci_chi_extra = 0.5 * vars.h[i][j] * (0.5 * h_UU_ww *
        dhww_dot_dchi + hUU_dchi_cartoon);

            out.LL[i][j] = ricci.LL[i][j] + nS * (ricci_chi_extra + vars.chi *
        ricci_tilde_extra) / vars.chi;

        }*/

        // Components of Ricci along the extra dimensions
        data_t dhww_cartoon = one_over_cartoon_coord * d1.hww[dI];
        data_t Gam_cartoon = one_over_cartoon_coord * vars.Gamma[dI];
        data_t reg12 = one_over_cartoon_coord2 * (h_UU[dI][dI] * vars.hww - 1.);

        data_t ricci_tilde_ww = 0;
        ricci_tilde_ww =
            -0.5 * nS * h_UU_ww * dhww_cartoon + vars.hww * Gam_cartoon + reg12;
        FOR(i)
        {
            // ricci_tilde_ww += 0.5 * (vars.Gamma[i] - 2 * Z_over_chi[i]) *
            // d1.hww[i];
            ricci_tilde_ww += 0.5 * chris_contracted[i] * d1.hww[i];
            FOR(j)
            {
                ricci_tilde_ww +=
                    0.5 * h_UU[i][j] *
                    (-d2.hww[i][j] + h_UU_ww * d1.hww[i] * d1.hww[j]);
            }
        }

        Tensor<2, data_t> covdtilde2chi;
        FOR(k, l)
        {
            covdtilde2chi[k][l] = d2.chi[k][l];
            FOR(m) { covdtilde2chi[k][l] -= chris.ULL[m][k][l] * d1.chi[m]; }
        }

        data_t boxtildechi = 0;
        data_t trdhhwwdch = 0;
        data_t tr_dch_dch = 0;
        FOR(k, l)
        {
            boxtildechi += h_UU[k][l] * covdtilde2chi[k][l];
            tr_dch_dch += h_UU[k][l] * d1.chi[k] * d1.chi[l];
            trdhhwwdch += h_UU[k][l] * d1.hww[k] * d1.chi[l];
        }

        // Debugging starts here. Alternative computation of Ricci
        Tensor<2, data_t> c_ri_chi;
        FOR(i, j)
        {
            c_ri_chi[i][j] =
                0.5 * vars.h[i][j] *
                    (boxtildechi +
                     nS * (0.5 * h_UU_ww * trdhhwwdch + hUU_dchi_cartoon) -
                     0.5 * GR_SPACEDIM * tr_dch_dch / vars.chi) +
                0.5 * (GR_SPACEDIM - 2) *
                    (covdtilde2chi[i][j] -
                     0.5 * d1.chi[i] * d1.chi[j] / vars.chi);
        }

        Tensor<2, data_t> c_ri;
        FOR(i, j)
        {
            c_ri[i][j] = c_ri_chi[i][j] + vars.chi * c_ri_hh[i][j];

            data_t z_terms = 0;
            FOR(k)
            {
                z_terms +=
                    Z_over_chi[k] *
                    (vars.h[i][k] * d1.chi[j] + vars.h[j][k] * d1.chi[i] -
                     vars.h[i][j] * d1.chi[k] + d1.h[i][j][k] * vars.chi);
            }
            c_ri[i][j] += z_terms;
            out.LL[i][j] = c_ri[i][j] / vars.chi;
        }
        // Debugging ends here

        data_t ricci_chi_ww = 0;
        ricci_chi_ww = 0.5 * vars.hww *
                       (boxtildechi +
                        (2 * GR_SPACEDIM - d - 2) *
                            (0.5 * h_UU_ww * trdhhwwdch + hUU_dchi_cartoon) -
                        0.5 * GR_SPACEDIM * tr_dch_dch / vars.chi);

        // out.LLww = (ricci_chi_ww + vars.chi * ricci_tilde_ww) / vars.chi;

        // Debugging: add Z terms to ww components
        data_t c_ri_ww = ricci_chi_ww + vars.chi * ricci_tilde_ww;
        FOR(i)
        {
            c_ri_ww -=
                Z_over_chi[i] * (vars.hww * d1.chi[i] - vars.chi * d1.hww[i]);
            out.LLww = c_ri_ww / vars.chi;
        }

        // Re-compute Ricci scalar
        out.scalar = vars.chi * (TensorAlgebra::compute_trace(out.LL, h_UU) +
                                 nS * h_UU_ww * out.LLww);

        // Debugging
        /*if(abs(x-m_center[0])<0.1 && abs(y-m_center[1])<0.1 &&
        abs((z-m_center[2])-5.)<m_dx)
        {
            std::cout.precision(15);
            //std::cout << " x = "<< x-m_center[0] << " y = "<< y-m_center[1] <<
        " z = "<< z-m_center[2] << " c_ri_xx = " << c_ri[0][0] << " c_ri_xy = "
        << c_ri[0][1] << " c_ri_xz = " << c_ri[0][2] << " c_ri_ww = " << c_ri_ww
        << " diff = " << abs(out.LL[0][0] - c_ri[0][0]/vars.chi) << std::endl;
            std::cout << " x = "<< x-m_center[0] << " y = "<< y-m_center[1] << "
        z = "<< z-m_center[2] << " Rxx * ch = " << out.LL[0][0] * vars.chi << "
        Rxy * ch = " << out.LL[0][1] * vars.chi << " Rxz * ch = " <<
        out.LL[0][2] * vars.chi << " Rww * ch = " << out.LLww * vars.chi <<
        std::endl;
        }*/

        return out;
    }

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    static ricciCartoon_t<data_t>
    compute_ricci(const vars_t<data_t> &vars,
                  const vars_t<Tensor<1, data_t>> &d1,
                  const diff2_vars_t<Tensor<2, data_t>> &d2,
                  const Tensor<2, data_t> &h_UU, const data_t &h_UU_ww,
                  const chris_t<data_t> &chris, const double &z)
    {
        Tensor<1, data_t> Z0 = 0.;
        return compute_ricci_Z(vars, d1, d2, h_UU, h_UU_ww, chris, Z0, z);
    }
};

#endif /* CCZ4CARTOONGEOMETRY_HPP_ */
