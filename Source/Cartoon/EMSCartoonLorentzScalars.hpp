/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMSCARTOONLORENTZSCALARS_HPP_
#define EMSCARTOONLORENTZSCALARS_HPP_

#include "ADMConformalVars.hpp" // needed for CCz4 and matter variables
#include "Cell.hpp"
// #include "EinsteinMaxwellDilatonField.hpp"
#include "MatterCCZ4.hpp" // need?
#include "CCZ4Vars.hpp" // need?
#include "FourthOrderDerivatives.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"
#include "EMSCouplingFunction.hpp"

//! Calculates the Noether Charge integrand values and the modulus of the
//! complex scalar field on the grid
template <class matter_t> class EMSCartoonLorentzScalars
{
    // Need matter variables and chi
    template <class data_t>
    using ADMVars = ADMConformalVars::VarsWithGauge<data_t>;
    template <class data_t>
    //using MatterVars = EinsteinMaxwellDilatonField<>::Vars<data_t>;
    using MatterVars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

    double m_dx;
    const std::array<double, CH_SPACEDIM> m_centre;
    CouplingFunction::params_t m_coupling_params;

  public:
    EMSCartoonLorentzScalars(double a_dx, std::array<double, CH_SPACEDIM> a_centre,
                                  CouplingFunction::params_t &a_coupling_params)
    : m_dx(a_dx), m_centre(a_centre), m_coupling_params(a_coupling_params) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const FourthOrderDerivatives m_deriv(m_dx);
        const auto adm_vars = current_cell.template load_vars<ADMVars>();
        const auto matter_vars = current_cell.template load_vars<MatterVars>();
        const auto adm_d1 = m_deriv.template diff1<ADMVars>(current_cell);
        const auto matter_d1 = m_deriv.template diff1<MatterVars>(current_cell);
        const auto matter_d2 = m_deriv.template diff2<MatterVars>(current_cell);

        // local variables for EMD coupling
        data_t alpha = m_coupling_params.alpha;
        data_t f0 = m_coupling_params.f0;
        data_t f1 = m_coupling_params.f1;
        data_t f2 = m_coupling_params.f2;

        // root minus g
        data_t root_minus_g = pow(adm_vars.chi, -1.5)*adm_vars.lapse;

        // inverse conformal metric
        const auto h_UU = TensorAlgebra::compute_inverse_sym(adm_vars.h);
        auto gamma_UU = h_UU;
        FOR2(i,j) gamma_UU[i][j] = h_UU[i][j]*adm_vars.chi;

        const auto chris = TensorAlgebra::compute_christoffel(adm_d1.h, h_UU);
        const auto chris_phys =
                           TensorAlgebra::compute_phys_chris(adm_d1.chi, adm_vars.chi,
                                                       adm_vars.h, h_UU, chris.ULL);

        // coordinates
        Coordinates<data_t> coords(current_cell, m_dx, m_centre);
        data_t x = coords.x, y = coords.y, z = coords.z;
        data_t r = sqrt(x*x + y*y + z*z);
        data_t safe_r = sqrt(r*r + 10e-20);
        data_t rho = sqrt(x*x + y*y);
        data_t safe_rho = sqrt(rho*rho + 10e-20);
        data_t costheta = z/safe_r;
        data_t sintheta = rho/safe_r;
        data_t cosphi = x/safe_rho;
        data_t sinphi = y/safe_rho;


        ////////////////////////////////////////////
        // calculate A_\mu A^\mu = -At^2 + a_i a^i
        // (At = A^\mu n_\mu = -\alpha A^t)
        ////////////////////////////////////////////
        data_t H2norm_damping = matter_vars.At * matter_vars.At
                              + matter_vars.Xi * matter_vars.Xi;


        ////////////////////////////////////////////
        // calculate F_{\mu\nu} F^{\mu\nu} = 2 B_i B^i - 2 E_i E^i
        ////////////////////////////////////////////

        data_t FF = 0.;
        Tensor<1, data_t, 3> Ei;
        Tensor<1, data_t, 3> Bi;
        data_t BB = 0.; // squared variables
        data_t EE = 0.; // squared variables

        Ei[0] = matter_vars.Ex;
        Ei[1] = matter_vars.Ey;
        Ei[2] = matter_vars.Ez;

        Bi[0] = matter_vars.ax;
        Bi[1] = matter_vars.ay;
        Bi[2] = matter_vars.az;


        FOR2(i, j)
        {
            EE += Ei[i] * Ei[j] * h_UU[i][j] * adm_vars.chi;
            BB += Bi[i] * Bi[j] * h_UU[i][j] * adm_vars.chi;
        }

        FF = 2. * BB - 2. * EE;


        ////////////////////////////////////////////
        // calculate hamiltonian of scalar field
        ////////////////////////////////////////////

        // data_t phi_hamiltonian = matter_vars.Pi * matter_vars.Pi;
        //
        // FOR1(i)
        // {
        //     phi_hamiltonian += 2. * matter_vars.Pi
        //                           * adm_vars.shift[i] * matter_d1.phi[i];
        // }
        //
        // FOR2(i,j)
        // {
        //     phi_hamiltonian += matter_d1.phi[i] * matter_d1.phi[j]
        //                      * h_UU[i][j] * adm_vars.chi;
        // }
        //
        // phi_hamiltonian *= root_minus_g;


        ////////////////////////////////////////////
        // Jacobeans
        ////////////////////////////////////////////

        // partial cartesian coords (i) by partial polar coords (i)
        // dxc_dxp[i][j]
        // i = {x,y,z}
        // j = {r,th,ph}
        Tensor<2, data_t, 3> dxc_dxp;
        dxc_dxp[0][0] = x/safe_r;
        dxc_dxp[1][0] = y/safe_r;
        dxc_dxp[2][0] = z/safe_r;
        dxc_dxp[0][1] = r * cosphi * costheta;
        dxc_dxp[1][1] = r * sinphi * costheta;
        dxc_dxp[2][1] = -r * sintheta;
        dxc_dxp[0][2] = -r * sinphi * sintheta;
        dxc_dxp[1][2] = r * cosphi * sintheta;
        dxc_dxp[2][2] = 0.; // dz/dphi=0

        // partial polar coords (i) by partial cartesian coords (i)
        // dxp_dxc[i][j]
        // i = {r,th,ph}
        // j = {x,y,z}
        Tensor<2, data_t, 3> dxp_dxc;
        dxp_dxc[0][0] = x/safe_r;
        dxp_dxc[0][1] = y/safe_r;
        dxp_dxc[0][2] = z/safe_r;
        dxp_dxc[1][0] = x * z / (safe_r * safe_r * safe_rho);
        dxp_dxc[1][1] = y * z / (safe_r * safe_r * safe_rho);
        dxp_dxc[1][2] = - rho / (safe_r * safe_r);
        dxp_dxc[2][0] = - y / (safe_rho * safe_rho);
        dxp_dxc[2][1] = x / (safe_rho * safe_rho);
        dxp_dxc[2][2] = 0.; // dphi/dz=0



        ////////////////////////////////////////////
        // maxwell constriants
        ////////////////////////////////////////////

        data_t H2norm_maxwell_constraints = 0.;
        Tensor<2, data_t, 3> DiBj;
        Tensor<2, data_t, 3> DiEj;
        data_t magnetic_constraint=0.;
        data_t electric_constraint=0.;
        data_t EDphi=0.; //E dot D phi

        // declare the default coupling values
        data_t f_of_phi = 0.0;
        data_t f_prime_of_phi = 0.0;
        data_t coupling_of_phi = 1.0;
        data_t local_phi = matter_vars.phi;

        compute_coupling_of_phi<data_t>(alpha,f0,f1,f2,f_of_phi,
                                f_prime_of_phi,coupling_of_phi,local_phi);

        FOR1(i)
        {
            DiBj[i][0] = matter_d1.ax[i];
            DiBj[i][1] = matter_d1.ay[i];
            DiBj[i][2] = matter_d1.az[i];

            DiEj[i][0] = matter_d1.Ex[i];
            DiEj[i][1] = matter_d1.Ey[i];
            DiEj[i][2] = matter_d1.Ez[i];
        }

        // partial derivs
        FOR2(i,j)
        {
            EDphi += Ei[i] * matter_d1.phi[j] * h_UU[i][j] * adm_vars.chi;
            electric_constraint += DiEj[i][j] * h_UU[i][j] * adm_vars.chi;
            magnetic_constraint += DiBj[i][j] * h_UU[i][j] * adm_vars.chi;
        }

        // covariant corrections
        FOR3(i,j,k)
        {
            magnetic_constraint += - chris_phys[k][i][j] * Bi[k]
                                     * h_UU[i][j] * adm_vars.chi;
            electric_constraint += - chris_phys[k][i][j] * Ei[k]
                                     * h_UU[i][j] * adm_vars.chi;
        }

        electric_constraint = coupling_of_phi * ( electric_constraint
                               - 2. * EDphi * f_prime_of_phi );

        H2norm_maxwell_constraints = electric_constraint * electric_constraint
                                    + magnetic_constraint * magnetic_constraint;



        ////////////////////////////////////////////
        // Mass and Charge Scalars
        ////////////////////////////////////////////

        // ADM mass calculation
        data_t g_rr = 0.;
        data_t ADM_scalar = 0.;
        Tensor<1, data_t, 3> dxdr; // equiv to sL, non normalised radial normal
        Tensor<1, data_t, 3> sU; // non-normalised unit upstairs radial vec
        dxdr[0] = x/safe_r;
        dxdr[1] = y/safe_r;
        dxdr[2] = z/safe_r;
        sU[0] = 0.;
        sU[1] = 0.;
        sU[2] = 0.;

        FOR2(i,j)
        {
            sU[i] += gamma_UU[i][j]*dxdr[j];
            g_rr += adm_vars.h[i][j]*dxdr[i]*dxdr[j]/adm_vars.chi;
        }
        FOR3(i,j,k)
        {
            ADM_scalar += sU[i]*gamma_UU[j][k]*(
                      (adm_d1.h[i][k][j]-adm_d1.h[j][k][i])/adm_vars.chi
                       -(adm_d1.chi[j]*adm_vars.h[i][k] -
                         adm_d1.chi[i]*adm_vars.h[j][k])*pow(adm_vars.chi,-2));                              ;
        }

        data_t Y_00 = 1./sqrt(4.*M_PI);

        // charge calculation
        data_t EUr = 0.; // E^r
        Tensor<1, data_t, 3> drdx;
        drdx[0] = x/safe_r;
        drdx[1] = y/safe_r;
        drdx[2] = z/safe_r;

        Tensor<2, data_t, 3> gamma_polar;
        FOR2(i,j)
        {
            EUr += Ei[i] * drdx[j] * gamma_UU[i][j];
            gamma_polar[i][j] = 0.;
        }

        FOR4(i,j,m,n)
        {
            gamma_polar[i][j] += adm_vars.h[m][n]/adm_vars.chi * dxc_dxp[m][i] * dxc_dxp[n][j];
        }
        data_t det_gamma_polar = TensorAlgebra::compute_determinant_sym(gamma_polar);
        data_t root_gamma_polar = sqrt(det_gamma_polar);

        data_t Q_scalar = EUr * coupling_of_phi * root_gamma_polar / sintheta;
        Q_scalar = Q_scalar/sqrt(2. * M_PI); // from weird def of Q in lagragean
        ADM_scalar = ADM_scalar * root_gamma_polar / (sintheta * 16. * M_PI);

        // store variables
        current_cell.store_vars(FF, c_mod_F);
        current_cell.store_vars(H2norm_maxwell_constraints, c_maxwell);
        // dividing by Y_00 makes the f_00 coeffcients actaully equal a regular
        // integral round teh circle
        current_cell.store_vars(ADM_scalar/Y_00, c_Mscalar);
        current_cell.store_vars(Q_scalar/Y_00, c_Qscalar);
    }
};

#endif /* EMSCARTOONLORENTZSCALARS_HPP_ */
