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
template <class coupling_t> class EMSCartoonLorentzScalars
{
    // Need matter variables and chi
    template <class data_t>
    using Vars = CCZ4CartoonVars::VarsWithGauge<data_t>;

    double m_dx;
    const std::array<double, CH_SPACEDIM> m_centre;
    CouplingFunction::params_t m_coupling_params;

  public:
    EMSCartoonLorentzScalars(double a_dx,
                            const std::array<double, CH_SPACEDIM> a_centre,
                                  CouplingFunction::params_t &a_coupling_params)
    : m_dx(a_dx), m_centre(a_centre), m_coupling_params(a_coupling_params) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const FourthOrderDerivatives m_deriv(m_dx);
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        //const auto d2 = m_deriv.template diff2<Vars>(current_cell);

        // local variables for EMD coupling
        data_t alpha = m_coupling_params.alpha;
        data_t f0 = m_coupling_params.f0;
        data_t f1 = m_coupling_params.f1;
        data_t f2 = m_coupling_params.f2;

        // pout() << "m_centre : " << m_centre[0] << ", " << m_centre[1] << std::endl;

        // root minus g
        data_t root_minus_g = pow(vars.chi, -1.5)*vars.lapse;

        // inverse conformal metric
        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        auto gamma_UU = h_UU;
        FOR2(i,j) gamma_UU[i][j] = h_UU[i][j]*vars.chi;

        // coordinates
        Coordinates<data_t> coords(current_cell, m_dx, m_centre);
        data_t x = coords.x, y = coords.y, z = 0.0;
        data_t r = sqrt(x*x + y*y + z*z);
        data_t safe_r = sqrt(r*r + 10e-20);
        data_t rho = sqrt(x*x + y*y);
        data_t safe_rho = sqrt(rho*rho + 10e-20);
        data_t costheta = z/safe_r;
        data_t sintheta = rho/safe_r;
        data_t cosphi = x/safe_rho;
        data_t sinphi = y/safe_rho;
        data_t ooy = 1./(y+10e-20);


        ////////////////////////////////////////////
        // mod sqaured of both damping params
        ////////////////////////////////////////////
        data_t H2norm_damping = vars.Lambda * vars.Lambda
                              + vars.Xi * vars.Xi;


        ////////////////////////////////////////////
        // calculate F_{\mu\nu} F^{\mu\nu} = 2 B_i B^i - 2 E_i E^i
        ////////////////////////////////////////////

        data_t FF = 0.;
        Tensor<1, data_t, 3> Ei;
        Tensor<1, data_t, 3> Bi;
        data_t BB = 0.; // squared variables
        data_t EE = 0.; // squared variables

        Ei[0] = vars.Ex;
        Ei[1] = vars.Ey;
        Ei[2] = vars.Ez;

        Bi[0] = vars.Bx;
        Bi[1] = vars.By;
        Bi[2] = vars.Bz;


        FOR2(i, j)
        {
            EE += Ei[i] * Ei[j] * h_UU[i][j] * vars.chi;
            BB += Bi[i] * Bi[j] * h_UU[i][j] * vars.chi;
        }
        EE += Ei[2]*Ei[2]*vars.chi/vars.hww;
        BB += Bi[2]*Bi[2]*vars.chi/vars.hww;

        FF = 2. * BB - 2. * EE;



        ////////////////////////////////////////////
        // scalar radiation
        ////////////////////////////////////////////
        data_t scalar_radiation = 0.;
        Tensor<1, data_t, 3> p_phi, p_phi_U;
        FOR(i)
        {
            // co-variant
            p_phi[i] = -2. * vars.Pi * d1.phi[i];
            p_phi_U[i] = 0.;
        }
        p_phi_U[2]=0;
        p_phi[2]=0.;

        FOR2(i,j)
        {
            //contra-variant
            p_phi_U[i] += p_phi[j] * gamma_UU[i][j];
        }
        // z component is zero anyway
        p_phi_U[2] = p_phi[2] * vars.hww * vars.chi;

        // radial radiation component
        // this line implicitly multiplied by r for the volume element
        scalar_radiation = (p_phi_U[0] * x + p_phi_U[1] * y) * pow(vars.chi,-1.5);


        ////////////////////////////////////////////
        // calculate hamiltonian of scalar field
        ////////////////////////////////////////////

        // data_t phi_hamiltonian = vars.Pi * vars.Pi;
        //
        // FOR1(i)
        // {
        //     phi_hamiltonian += 2. * vars.Pi
        //                           * vars.shift[i] * d1.phi[i];
        // }
        //
        // FOR2(i,j)
        // {
        //     phi_hamiltonian += d1.phi[i] * d1.phi[j]
        //                      * h_UU[i][j] * vars.chi;
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
        Tensor<1, data_t, 3> drootgammathing;
        FOR1(i) drootgammathing[i] = -1.5*d1.chi[i]/vars.chi;
        drootgammathing[2] = 0.;

        Tensor<3, data_t, 3> dgamma3d = {0.};
        FOR3(a,b,c)
        {
            dgamma3d[a][b][c] = (d1.h[b][c][a]
                                - d1.chi[a]*vars.h[b][c]/vars.chi)/vars.chi;
        }
        dgamma3d[0][2][2] = (d1.hww[0]
                            - d1.chi[0]*vars.hww/vars.chi)/vars.chi;
        dgamma3d[1][2][2] = (d1.hww[1]
                            - d1.chi[1]*vars.hww/vars.chi)/vars.chi;
        dgamma3d[2][0][2] = ooy * vars.h[0][1] / vars.chi;
        dgamma3d[2][2][0] = dgamma3d[2][0][2];
        dgamma3d[2][1][2] = ooy * (vars.h[1][1] - vars.hww)/vars.chi;
        dgamma3d[2][2][1] = dgamma3d[2][1][2];

        Tensor<2, data_t, 3> gammaUU_3d = {0.};
        FOR2(i,j) gammaUU_3d[i][j] = vars.chi * h_UU[i][j];
        gammaUU_3d[2][2] = vars.chi / vars.hww;

        Tensor<1, data_t, 3> div_gamma_inv = {0.}; // parial_a gamma^{ad}
        for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
        for (int c = 0; c < 3; c++) {
        for (int d = 0; d < 3; d++) {
            div_gamma_inv[d] += - gammaUU_3d[a][b] * gammaUU_3d[c][d]
                                * dgamma3d[a][b][c];
        }}}}

        // declare the default coupling values
        data_t f_of_phi = 0.0;
        data_t f_prime_of_phi = 0.0;
        data_t coupling_of_phi = 1.0;
        data_t local_phi = vars.phi;

        compute_coupling_of_phi<data_t>(alpha,f0,f1,f2,f_of_phi,
                                f_prime_of_phi,coupling_of_phi,local_phi);

        FOR1(i)
        {
            DiBj[i][0] = d1.Bx[i];
            DiBj[i][1] = d1.By[i];
            DiBj[i][2] = d1.Bz[i];

            DiEj[i][0] = d1.Ex[i];
            DiEj[i][1] = d1.Ey[i];
            DiEj[i][2] = d1.Ez[i];
        }
        DiBj[2][0] = 0.;
        DiBj[2][1] = - ooy * vars.Bz;
        DiBj[2][2] = ooy * vars.By;

        DiEj[2][0] = 0.;
        DiEj[2][1] = - ooy * vars.Ez;
        DiEj[2][2] = ooy * vars.Ey;

        // partial derivs
        FOR2(i,j)
        {
            EDphi += Ei[i] * d1.phi[j] * h_UU[i][j] * vars.chi;
            electric_constraint += DiEj[i][j] * h_UU[i][j] * vars.chi;
            magnetic_constraint += DiBj[i][j] * h_UU[i][j] * vars.chi;
            electric_constraint += -1.5 * h_UU[i][j] * Ei[i] * d1.chi[j];
            magnetic_constraint += -1.5 * h_UU[i][j] * Bi[i] * d1.chi[j];
        }
        FOR1(i)
        {
            electric_constraint += Ei[i] * div_gamma_inv[i];
            magnetic_constraint += Bi[i] * div_gamma_inv[i];
        }
        // cartoon terms
        electric_constraint += ooy*vars.Ey*vars.chi/vars.hww;
        magnetic_constraint += ooy*vars.By*vars.chi/vars.hww;


        electric_constraint = coupling_of_phi * ( electric_constraint
                               - 2. * EDphi * f_prime_of_phi );

        H2norm_maxwell_constraints = electric_constraint * electric_constraint
                                    + magnetic_constraint * magnetic_constraint;


        // electric_constraint = coupling_of_phi * ( electric_constraint
        //                        - 2. * EDphi * f_prime_of_phi );
        //
        // H2norm_maxwell_constraints = electric_constraint;

        ////////////////////////////////////////////
        // Mass and Charge Scalars
        ////////////////////////////////////////////

        // ADM mass calculation
        // data_t g_rr = 0.;
        data_t ADM_scalar = 0.;
        // Tensor<1, data_t, 3> dxdr; // equiv to sL, non normalised radial normal
        // Tensor<1, data_t, 3> sU; // non-normalised unit upstairs radial vec
        // dxdr[0] = x/safe_r;
        // dxdr[1] = y/safe_r;
        // dxdr[2] = z/safe_r;
        // sU[0] = 0.;
        // sU[1] = 0.;
        // sU[2] = 0.;
        //
        // for (int i = 0; i < 3; i++) {
        // for (int j = 0; j < 3; j++) {
        //     sU[i] += gamma_UU[i][j]*dxdr[j];
        //     g_rr += vars.h[i][j]*dxdr[i]*dxdr[j]/vars.chi;
        // }}
        // for (int i = 0; i < 3; i++) {
        // for (int j = 0; j < 3; j++) {
        // for (int k = 0; k < 3; k++) {
        //     ADM_scalar += sU[i]*gamma_UU[j][k]*(
        //               (d1.h[i][k][j]-d1.h[j][k][i])/vars.chi
        //                -(d1.chi[j]*vars.h[i][k] -
        //                  d1.chi[i]*vars.h[j][k])*pow(vars.chi,-2));                              ;
        // }}}

        data_t Y_00 = 1./sqrt(4.*M_PI);
        data_t four_pi = M_PI*4.;

        // charge calculation
        data_t EUr = 0.; // E^r
        data_t drdchi = 0.; // partial_r chi
        Tensor<1, data_t, 3> drdx;
        drdx[0] = x/safe_r;
        drdx[1] = y/safe_r;
        drdx[2] = z/safe_r;

        FOR2(i,j)
        {
            EUr += Ei[i] * drdx[j] * vars.chi * h_UU[i][j];
        }
        FOR(i)
        {
            drdchi += drdx[i] * d1.chi[i];
        }
        //cartoon terms = 0


        data_t Q_scalar = four_pi * coupling_of_phi * r * r * pow(vars.chi,-1.5) * EUr;
        Q_scalar = Q_scalar / sqrt(2. * M_PI); // from definition of E in Lagrangean

        // useful test line to see where centre is
        // Q_scalar = x*x + y*y;

        // ADM_scalar = -r^2 d_psi/d_r from robins thesis
        ADM_scalar = 0.5*r*r*pow(vars.chi,-1.5)*drdchi;

        // store variables
        current_cell.store_vars(FF, c_mod_F);
        current_cell.store_vars(scalar_radiation, c_phi_rad);
        // dividing by Y_00 makes the f_00 coeffcients actaully equal a regular
        // integral round teh circle

        // i think the 4pi here is a complication of the 2d code
        // it was calibrated by integrating 1 over a semi-circle centred y=0, x=L/2
        current_cell.store_vars(ADM_scalar/Y_00/four_pi, c_Mscalar);
        current_cell.store_vars(Q_scalar/Y_00/four_pi, c_Qscalar);
    }
};

#endif /* EMSCARTOONLORENTZSCALARS_HPP_ */
