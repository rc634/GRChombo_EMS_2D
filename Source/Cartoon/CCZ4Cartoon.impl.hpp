/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CCZ4CARTOON_HPP_)
#error "This file should only be included through CCZ4Cartoon.hpp"
#endif

#ifndef CCZ4CARTOON_IMPL_HPP_
#define CCZ4CARTOON_IMPL_HPP_

#define COVARIANTZ4
#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"
#include <iostream>
#include <stdio.h>

template <class gauge_t, class deriv_t, class coupling_t>
inline CCZ4Cartoon<gauge_t, deriv_t, coupling_t>::CCZ4Cartoon(
        params_t params, double dx, double sigma,
        coupling_t a_coupling, double a_G_Newton, int formulation,
        double cosmological_constant)
        : CCZ4RHS<gauge_t, deriv_t>(params, dx, sigma, formulation,
                                    cosmological_constant),
          m_coupling(a_coupling), m_G_Newton(a_G_Newton)
{
}

template <class data_t, template <typename> class vars_t>
emtensorCartoon_t<data_t>
compute_EMS_EM_tensor(const vars_t<data_t> &vars,
                     const vars_t<Tensor<1, data_t>> &d1,
const Tensor<2, data_t> &h_UU, const data_t &h_UU_ww,
const chris_t<data_t> &chris, const int &nS, const CouplingFunction &a_coupling)
{
emtensorCartoon_t<data_t> out;

// set coupling functions
data_t coupling =0., f_coupling=0., fprime=0.;
a_coupling.compute_coupling(f_coupling, fprime, coupling, vars);

// 3D B and E fields
Tensor<1, data_t, 3> E, B;
E[0] = vars.Ex;
E[1] = vars.Ey;
E[2] = vars.Ez;
B[0] = vars.Bx;
B[1] = vars.By;
B[2] = vars.Bz;

// epsilon tensor
Tensor<3, data_t, 3> epsLLL = {0.};
data_t eps012 = pow(vars.chi, -1.5);
epsLLL[0][1][2] = eps012;
epsLLL[0][2][1] = -eps012;
epsLLL[1][2][0] = eps012;
epsLLL[1][0][2] = -eps012;
epsLLL[2][0][1] = eps012;
epsLLL[2][1][0] = -eps012;

// 3d spatial inverse metric
Tensor<2, data_t, 3> gammaUU_3d = {0.};
FOR2(i,j)
{
    gammaUU_3d[i][j] = h_UU[i][j] * vars.chi;
}
gammaUU_3d[2][2] = h_UU_ww * vars.chi;

// eps_{ijk} E^j B^k
Tensor<1, data_t, 3> epsEB = {0.};
for (int i=0; i<3; i++) {
for (int j=0; j<3; j++) {
for (int k=0; k<3; k++) {
for (int m=0; m<3; m++) {
for (int n=0; n<3; n++) {
    epsEB[i] += epsLLL[i][j][k] * gammaUU_3d[j][m] * E[m]
                                * gammaUU_3d[k][n] * B[n];
}}}}}

// some bits we need
data_t Dphi_sqr = 0.;
data_t BB = 0.;
data_t EE = 0.;
FOR2(i,j)
{
    Dphi_sqr += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j];
    BB += vars.chi * h_UU[i][j] * B[i] * B[j];
    EE += vars.chi * h_UU[i][j] * E[i] * E[j];
}
// cartoon terms
BB += vars.chi * h_UU_ww * B[2] * B[2];
EE += vars.chi * h_UU_ww * E[2] * E[2];

// stress tensor components
out.rho = vars.Pi * vars.Pi + Dphi_sqr + coupling * (BB + EE);

FOR(i)
{
    out.Si[i] = 2. * (vars.Pi * d1.phi[i] + coupling * epsEB[i]);
}

FOR(i, j)
{
    out.Sij[i][j] = 2. * d1.phi[i] * d1.phi[j]
                       - vars.h[i][j] * (Dphi_sqr - vars.Pi * vars.Pi)/vars.chi
                    + coupling * (-2. * (E[i] * E[j] + B[i] * B[j])
                                      + vars.h[i][j] * (BB + EE)/vars.chi );
}



// Sww
out.Sww = - (vars.hww/vars.chi) * (Dphi_sqr - vars.Pi * vars.Pi)
                + coupling * (-2. * (vars.Ez*vars.Ez + vars.Bz*vars.Bz)
                                  + (vars.hww/vars.chi) * (BB + EE) );

// S
out.S = TensorAlgebra::compute_trace(out.Sij, h_UU);
//cartoon term
out.S += h_UU_ww * out.Sww;
// mult by chi as gamma_UU = chi * h_UU, used in trace
out.S *= vars.chi;

return out;
}

template <class gauge_t, class deriv_t, class coupling_t>
template <class data_t>
void CCZ4Cartoon<gauge_t, deriv_t, coupling_t>::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
            this->m_deriv.template advection<Vars>(current_cell, vars.shift);

    Coordinates<data_t> coords(current_cell, this->m_deriv.m_dx);

    Vars<data_t> rhs;
    CCZ4Cartoon<gauge_t, deriv_t>::rhs_equation(rhs, vars, d1, d2, advec,
                                                coords.y);

    this->m_deriv.add_dissipation(rhs, current_cell, this->m_sigma);

    current_cell.store_vars(rhs); // Write the rhs into the output FArrayBox
}

template <class gauge_t, class deriv_t, class coupling_t>
template <class data_t, template <typename> class vars_t,
        template <typename> class diff2_vars_t>
void CCZ4Cartoon<gauge_t, deriv_t, coupling_t>::rhs_equation(
        vars_t<data_t> &rhs, const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1,
const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
const double &cartoon_coord) const
{
using namespace TensorAlgebra;
const int dI = CH_SPACEDIM - 1;
const int nS =
        GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
const double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
const double two_over_gr_spacedim = 2. * one_over_gr_spacedim;

auto h_UU = compute_inverse_sym(vars.h);
auto h_UU_ww = 1. / vars.hww;
auto chris = compute_christoffel(d1.h, h_UU);
const double one_over_cartoon_coord = 1. / cartoon_coord;
const double one_over_cartoon_coord2 = 1. / (cartoon_coord * cartoon_coord);

Tensor<1, data_t> chris_ww;
FOR(i)
        {
                chris_ww[i] =
                        one_over_cartoon_coord * (delta(i, dI) - h_UU[i][dI] * vars.hww);
        FOR(j) chris_ww[i] -= 0.5 * h_UU[i][j] * d1.hww[j];
        }
Tensor<1, data_t> chris_contracted; //!< includes the higher D contributions
FOR(i)
        {
                chris_contracted[i] = chris.contracted[i] + nS * h_UU_ww * chris_ww[i];
        }

Tensor<1, data_t> Z_over_chi;
Tensor<1, data_t> Z;
if (this->m_formulation == CCZ4::USE_BSSN)
{
FOR(i) Z_over_chi[i] = 0.0;
}
else
{
FOR(i) Z_over_chi[i] = 0.5 * (vars.Gamma[i] - chris_contracted[i]);
}
FOR(i) Z[i] = vars.chi * Z_over_chi[i];

auto ricci = CCZ4CartoonGeometry::compute_ricci_Z(
        vars, d1, d2, h_UU, h_UU_ww, chris, Z_over_chi, cartoon_coord);

data_t divshift = compute_trace(d1.shift) +
                  nS * vars.shift[dI] /
                  cartoon_coord; //!< includes the higher D contribution
data_t divshift_w = compute_trace(d1.shift) -
                    CH_SPACEDIM * vars.shift[dI] /
                    cartoon_coord; //!< combination that appears in the
//!< rhs of the ww components
data_t Z_dot_d1lapse = compute_dot_product(Z, d1.lapse);
data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);
data_t dlapse_dot_dhww = compute_dot_product(d1.lapse, d1.hww, h_UU);

data_t hUU_dlapse_cartoon = 0;
FOR(i)
        {
                hUU_dlapse_cartoon +=
                        one_over_cartoon_coord * h_UU[dI][i] * d1.lapse[i];
        }

Tensor<1, data_t> reg_03, reg_04;
FOR(i)
        {
                reg_03[i] = one_over_cartoon_coord * d1.shift[i][dI] -
                            one_over_cartoon_coord2 * delta(i, dI) * vars.shift[dI];
        reg_04[i] = one_over_cartoon_coord * d1.shift[dI][i] -
        one_over_cartoon_coord2 * delta(i, dI) * vars.shift[dI];
        }

Tensor<2, data_t> covdtilde2lapse;
Tensor<2, data_t> covd2lapse; // NOTE: we compute chi * D_i D_j lapse
FOR(k, l)
{
covdtilde2lapse[k][l] = d2.lapse[k][l];
FOR(m) { covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m]; }
covd2lapse[k][l] =
vars.chi * covdtilde2lapse[k][l] +
0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
vars.h[k][l] * dlapse_dot_dchi);
}
data_t covd2lapse_ww =
        vars.chi * (0.5 * dlapse_dot_dhww + vars.hww * hUU_dlapse_cartoon) -
        0.5 * vars.hww * dlapse_dot_dchi;

/*data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
FOR(i)
{
    tr_covd2lapse -= vars.chi * chris.contracted[i] * d1.lapse[i];
    FOR(j)
    {
        tr_covd2lapse += h_UU[i][j] * (vars.chi * d2.lapse[i][j] +
                                       d1.lapse[i] * d1.chi[j]);
    }
}*/
data_t tr_covd2lapse = compute_trace(covd2lapse, h_UU);
tr_covd2lapse += nS * h_UU_ww * covd2lapse_ww;

// Make A_{ij} trace free
/*data_t trA = nS * h_UU_ww * vars.Aww;
FOR(i,j)
{
    trA += h_UU[i][j] * vars.A[i][j];
}*/
data_t trA = compute_trace(vars.A, h_UU);
trA += nS * h_UU_ww * vars.Aww;

Tensor<2, data_t> A_TF;
FOR(i, j)
{
A_TF[i][j] = vars.A[i][j] - trA * vars.h[i][j] * one_over_gr_spacedim;
}
data_t Aww_TF;
Aww_TF = vars.Aww - trA * vars.hww * one_over_gr_spacedim;

Tensor<2, data_t> A_UU = raise_all(A_TF, h_UU);
data_t A_UU_ww = h_UU_ww * h_UU_ww * Aww_TF;
// A^{ij} A_{ij}. - Note the abuse of the compute trace function.
data_t tr_A2 = compute_trace(A_TF, A_UU) + nS * A_UU_ww * Aww_TF;

rhs.chi = advec.chi + two_over_gr_spacedim * vars.chi *
                      (vars.lapse * vars.K - divshift);

FOR(i, j)
{
rhs.h[i][j] = advec.h[i][j] - 2.0 * vars.lapse * A_TF[i][j] -
two_over_gr_spacedim * vars.h[i][j] * divshift;
FOR(k)
        {
                rhs.h[i][j] +=
                        vars.h[k][i] * d1.shift[k][j] + vars.h[k][j] * d1.shift[k][i];
        }
}

rhs.hww = advec.hww - 2.0 * vars.lapse * Aww_TF -
          two_over_gr_spacedim * vars.hww * divshift_w;

Tensor<2, data_t> Adot_TF;
FOR(i, j)
{
Adot_TF[i][j] =
-covd2lapse[i][j] + vars.chi * vars.lapse * ricci.LL[i][j];
}
data_t Adot_TF_ww = -covd2lapse_ww + vars.chi * vars.lapse * ricci.LLww;

data_t trAdot = compute_trace(Adot_TF, h_UU) + nS * h_UU_ww * Adot_TF_ww;
// make_trace_free(Adot_TF, vars.h, h_UU);
FOR(i, j) { Adot_TF[i][j] -= one_over_gr_spacedim * vars.h[i][j] * trAdot; }
Adot_TF_ww -= one_over_gr_spacedim * vars.hww * trAdot;

FOR(i, j)
{
rhs.A[i][j] = advec.A[i][j] + Adot_TF[i][j] +
vars.A[i][j] * (vars.lapse * (vars.K - 2 * vars.Theta) -
two_over_gr_spacedim * divshift);
FOR(k)
        {
                rhs.A[i][j] +=
                        A_TF[k][i] * d1.shift[k][j] + A_TF[k][j] * d1.shift[k][i];
        FOR(l)
        {
            rhs.A[i][j] -=
                    2 * vars.lapse * h_UU[k][l] * A_TF[i][k] * A_TF[l][j];
        }
        }
}

rhs.Aww = advec.Aww + Adot_TF_ww +
          Aww_TF * (vars.lapse *
                    ((vars.K - 2 * vars.Theta) - 2 * h_UU_ww * Aww_TF) -
                    two_over_gr_spacedim * divshift_w);

#ifdef COVARIANTZ4
data_t kappa1_lapse = this->m_params.kappa1;
#else
data_t kappa1_lapse = this->m_params.kappa1 * vars.lapse;
#endif

rhs.Theta =
advec.Theta +
0.5 * vars.lapse *
(ricci.scalar - tr_A2 +
 ((GR_SPACEDIM - 1.0) * one_over_gr_spacedim) * vars.K * vars.K -
 2 * vars.Theta * vars.K) -
0.5 * vars.Theta * kappa1_lapse *
((GR_SPACEDIM + 1) + this->m_params.kappa2 * (GR_SPACEDIM - 1)) -
Z_dot_d1lapse;

rhs.Theta += -vars.lapse * this->m_cosmological_constant;

rhs.K =
advec.K +
vars.lapse * (ricci.scalar + vars.K * (vars.K - 2 * vars.Theta)) -
kappa1_lapse * GR_SPACEDIM * (1 + this->m_params.kappa2) * vars.Theta -
tr_covd2lapse;
rhs.K += -2 * vars.lapse * GR_SPACEDIM / ((double)GR_SPACEDIM - 1.) *
this->m_cosmological_constant;

Tensor<1, data_t> Gammadot;
FOR(i)
        {
                Gammadot[i] =
                        two_over_gr_spacedim *
                        (divshift * (chris_contracted[i] +
                                     2 * this->m_params.kappa3 * Z_over_chi[i]) -
                         2 * vars.lapse * vars.K * Z_over_chi[i]) -
                        2 * kappa1_lapse * Z_over_chi[i] +
                        nS *
                        (2. * vars.lapse * chris_ww[i] * A_UU_ww + h_UU_ww * reg_03[i]);
        FOR(j)
        {
            Gammadot[i] +=
                    2 * h_UU[i][j] *
                    (vars.lapse * d1.Theta[j] - vars.Theta * d1.lapse[j]) -
                    2 * A_UU[i][j] * d1.lapse[j] -
                    vars.lapse * (two_over_gr_spacedim * (GR_SPACEDIM - 1.0) *
                                  h_UU[i][j] * d1.K[j] +
                                  GR_SPACEDIM * A_UU[i][j] * d1.chi[j] / vars.chi) -
                    (chris_contracted[j] +
                     2 * this->m_params.kappa3 * Z_over_chi[j]) *
                    d1.shift[i][j] +
                    (nS * (GR_SPACEDIM - 2.0) * one_over_gr_spacedim) * h_UU[i][j] *
                    reg_04[j];
            FOR(k)
            {
                Gammadot[i] +=
                        2 * vars.lapse * chris.ULL[i][j][k] * A_UU[j][k] +
                        h_UU[j][k] * d2.shift[i][j][k] +
                        ((GR_SPACEDIM - 2.0) * one_over_gr_spacedim) * h_UU[i][j] *
                        d2.shift[k][j][k];
            }
        }
        rhs.Gamma[i] = advec.Gamma[i] + Gammadot[i];
        }

this->m_gauge.rhs_gauge(rhs, vars, d1, d2, advec);

////////////////////////
// Stress Tensor Calc
///////////////////////

// calculate coupling
data_t coupling =0., f_coupling=0., fprime=0.;

m_coupling.compute_coupling(f_coupling, fprime, coupling, vars);

emtensorCartoon_t<data_t> emtensor =
        compute_EMS_EM_tensor(vars, d1, h_UU, h_UU_ww, chris, nS, m_coupling);

Tensor<2, data_t> Sij_TF;
data_t Sww_TF;

FOR(i, j)
{
    Sij_TF[i][j] = emtensor.Sij[i][j]
                  - one_over_gr_spacedim * emtensor.S * vars.h[i][j] / vars.chi;
}
Sww_TF = emtensor.Sww - one_over_gr_spacedim * emtensor.S * vars.hww / vars.chi;

FOR(i, j)
{
rhs.A[i][j] +=
-8.0 * M_PI * m_G_Newton * vars.chi * vars.lapse * Sij_TF[i][j];
}

rhs.Aww += -8.0 * M_PI * m_G_Newton * vars.chi * vars.lapse * Sww_TF;

rhs.K +=
4.0 * M_PI * m_G_Newton * vars.lapse * (emtensor.S - 3 * emtensor.rho);
rhs.Theta += -8.0 * M_PI * m_G_Newton * vars.lapse * emtensor.rho;

FOR(i)
{
    data_t matter_term_Gamma = 0.0;
    FOR(j)
    {
        matter_term_Gamma += -16.0 * M_PI * m_G_Newton * vars.lapse *
                             h_UU[i][j] * emtensor.Si[j];
    }

    rhs.Gamma[i] += matter_term_Gamma;
    rhs.B[i] += matter_term_Gamma;
}

////////////////////////////////////////
// Matter evolution equations
///////////////////////////////////////

// some conveniences and damping parameters
data_t ooy = sqrt(1./(pow(1./one_over_cartoon_coord,2)+10e-20)); // 1/y
const int iy = dI; // index of y coord
data_t kappa_B = 1.0, kappa_E = 1.0;

Tensor<1, data_t, 3> E;
Tensor<1, data_t, 3> B;
E[0] = vars.Ex;
E[1] = vars.Ey;
E[2] = vars.Ez;
B[0] = vars.Bx;
B[1] = vars.By;
B[2] = vars.Bz;

Tensor<3, data_t, 3> epsLLL = {0.};
data_t eps012 = pow(vars.chi, -1.5);
epsLLL[0][1][2] = eps012;
epsLLL[0][2][1] = -eps012;
epsLLL[1][2][0] = eps012;
epsLLL[1][0][2] = -eps012;
epsLLL[2][0][1] = eps012;
epsLLL[2][1][0] = -eps012;

// some derivatives of tensors
Tensor<3, data_t, 3> dgammainv = {0.}; // d[i] gammainv ^[j,k]
Tensor<3, data_t, 3> dgamma3d = {0.}; // d[i] gamma [j,k]
Tensor<2, data_t, 3> dE; // d[i] E[j] = dE[i][j]
Tensor<2, data_t, 3> dB; // d[i] B[j] = dB[i][j]
Tensor<2, data_t> Kij; // K_{ij}
data_t Kzz = vars.Aww / vars.chi + vars.K * vars.hww / 3.;
FOR2(i,j)
{
    Kij[i][j] = vars.A[i][j]/vars.chi + vars.K*vars.h[i][j]/3.;
}

FOR(i)
{
    dE[i][0] = d1.Ex[i];
    dE[i][1] = d1.Ey[i];
    dE[i][2] = d1.Ez[i];
    dB[i][0] = d1.Bx[i];
    dB[i][1] = d1.By[i];
    dB[i][2] = d1.Bz[i];
}
//cartoon terms for z derivatives
dE[2][0] = 0.;
dE[2][1] = - ooy * vars.Ez;
dE[2][2] = ooy * vars.Ey;
dB[2][0] = 0.;
dB[2][1] = - ooy * vars.Bz;
dB[2][2] = ooy * vars.By;



// compute the terms that are the cross-product terms in maxwell
// i.e. epsilon tensor terms,
// note E(B) term containts B(E) fields (opposite)
Tensor<1, data_t, 3> d_lapse_3d;
Tensor<1, data_t, 3> d_phi_3d;
FOR(i)
{
    d_lapse_3d[i] = d1.lapse[i];
    d_phi_3d[i] = d1.phi[i];
}
d_lapse_3d[2] = 0.;
d_phi_3d[2] = 0.;

Tensor<2, data_t, 3> h_UU_3d = {0.};
FOR2(i,j)
{
    h_UU_3d = h_UU[i][j];
}
h_UU_3d[2][2] = h_UU_ww;


FOR3(i,j,k)
{
    dgamma3d[i][j][k] = d1.h[j][k][i] / vars.chi
                       - vars.h[j][k] * d1.chi[i] * pow(vars.chi , -2);
}
dgamma3d[0][2][2] = d1.hww[0] / vars.chi
                       - vars.hww * d1.chi[0] * pow(vars.chi , -2);
dgamma3d[1][2][2] = d1.hww[1] / vars.chi
                      - vars.hww * d1.chi[1] * pow(vars.chi , -2);

dgamma3d[2][0][2] = ooy * vars.h[0][1] / vars.chi;
dgamma3d[2][2][0] = dgamma3d[2][0][2];
dgamma3d[2][1][2] = ooy * (vars.h[1][1] - vars.hww)/vars.chi;
dgamma3d[2][2][1] = dgamma3d[2][1][2];

for (int i=0; i<3; i++){
for (int a=0; a<3; a++){
for (int b=0; b<3; b++){
for (int c=0; c<3; c++){
for (int d=0; d<3; d++){
    dgammainv[i][a][d] += - vars.chi * vars.chi
                              * h_UU_3d[c][d] * h_UU_3d[a][b]
                              * dgamma3d[i][b][c];
}}}}}


Tensor<2, data_t, 3> the_E_term;
Tensor<2, data_t, 3> the_B_term;
for (int k = 0; k <3; k++){
for (int j = 0; j <3; j++){
    the_E_term[k][j] = B[k] * d_lapse_3d[j] + vars.lapse * dB[j][k]
                            - 2. * fprime * B[k] * d_phi_3d[j];
    the_B_term[k][j] = - E[k] * d_lapse_3d[j] - vars.lapse * dE[j][k];
}}

// Make some objects we need to keep the equations readable
data_t FF = 0.;
data_t EdotDphi = 0.;
data_t divB3D = 0.;
data_t divE3D = 0.;
data_t EKx = 0., EKy = 0., EzKzz = vars.Ez * vars.chi * h_UU_ww * Kzz;
data_t BKx = 0., BKy = 0., BzKzz = vars.Bz * vars.chi * h_UU_ww * Kzz;
FOR2(i,j)
{

    EKx += E[i] * Kij[0][j] * vars.chi * h_UU[i][j];
    EKy += E[i] * Kij[1][j] * vars.chi * h_UU[i][j];
    BKx += B[i] * Kij[0][j] * vars.chi * h_UU[i][j];
    BKy += B[i] * Kij[1][j] * vars.chi * h_UU[i][j];

    EdotDphi += vars.chi * h_UU[i][j] * E[i] * d1.phi[j];

    divE3D += -3./2. * h_UU[i][j] * E[i] * d1.chi[j];
    divE3D += vars.chi * h_UU[i][j] * dE[i][j]; // the "dE term"

    divB3D += -3./2. * h_UU[i][j] * B[i] * d1.chi[j];
    divB3D += vars.chi * h_UU[i][j] * dB[i][j]; // the "dB term"

    FF += 2. * vars.chi * h_UU[i][j] * (B[i] * B[j] - E[i] * E[j]);
}

//cartoon terms
divE3D += ooy * h_UU_ww * vars.Ey * vars.chi; // the cartoon "dE term"
divB3D += ooy * h_UU_ww * vars.By * vars.chi; // the cartoon "dB term"
FF += 2. * vars.chi * h_UU_ww * (B[2] * B[2] - E[2] * E[2]);

// need 3d loop here as d_z gamma^zx =/=0 (even thought gamma^xz=0)
// same arguemtn for gamma^zy ...
for (int i = 0; i <3; i++){
for (int j = 0; j <3; j++){
  divE3D += dgammainv[i][i][j] * E[j];
  divB3D += dgammainv[i][i][j] * B[j];
}}




////////////////////////////////
// Lie Derivatives
////////////////////////////////

// Lie derivatives, uses d_i beta_j = d1.shift[j][i]
data_t LieBetaEx = advec.Ex
                    + vars.Ex * d1.shift[0][0] + vars.Ey * d1.shift[1][0];
data_t LieBetaEy = advec.Ey
                    + vars.Ex * d1.shift[0][1] + vars.Ey * d1.shift[1][1];
data_t LieBetaEz = advec.Ez + ooy * vars.Ez * vars.shift[iy];
data_t LieBetaBx = advec.Bx
                    + vars.Bx * d1.shift[0][0] + vars.By * d1.shift[1][0];
data_t LieBetaBy = advec.By
                    + vars.Bx * d1.shift[0][1] + vars.By * d1.shift[1][1];
data_t LieBetaBz = advec.Bz + ooy * vars.Bz * vars.shift[iy];


// The actual equations
rhs.phi = advec.phi - vars.lapse * vars.Pi;
rhs.Pi = advec.Pi + vars.lapse * vars.K * vars.Pi
                  - ooy * h_UU_ww * vars.chi * d1.phi[iy]
                  -0.5 * vars.lapse * fprime * coupling * FF;

rhs.Ex = LieBetaEx + vars.lapse * (vars.K * vars.Ex - 2. * EKx
                                   + d1.Xi[0]
                                   - 2. * fprime * vars.Pi * vars.Ex);
rhs.Ey = LieBetaEy + vars.lapse * (vars.K * vars.Ey - 2. * EKy
                                   + d1.Xi[1]
                                   - 2. * fprime * vars.Pi * vars.Ey);
rhs.Ez = LieBetaEz + vars.lapse * (vars.K * vars.Ez - 2. * EzKzz
                                   - 2. * fprime * vars.Pi * vars.Ez);
rhs.Bx = LieBetaBx + vars.lapse * (vars.K * vars.Bx - 2. * BKx
                                  + d1.Lambda[0]
                                  - 2. * fprime * vars.Pi * vars.Ex);
rhs.By = LieBetaBy + vars.lapse * (vars.K * vars.By - 2. * BKy
                                  + d1.Lambda[1]
                                  - 2. * fprime * vars.Pi * vars.Ey);
rhs.Bz = LieBetaBz + vars.lapse * (vars.K * vars.Bz - 2. * BzKzz);
rhs.Lambda = advec.Lambda + vars.lapse * ( divB3D - kappa_B * vars.Lambda);
rhs.Xi = advec.Xi + vars.lapse * coupling * (divE3D - 2. * fprime * EdotDphi)
                                          - vars.lapse * kappa_E * vars.Xi;

FOR (i)
{
    rhs.Pi += vars.chi * vars.lapse * vars.Gamma[i] * d1.phi[i];
}

// note the 0.5*d1.chi term comes from the conformal christoffel
FOR2  (i,j)
{
    rhs.Pi += h_UU[i][j]*( - vars.chi * vars.lapse * d2.phi[i][j]
                           + 0.5 * d1.chi[i] * vars.lapse * d1.phi[j]
                           - vars.chi * d1.lapse[i] * d1.phi[j]);
}

// explicit loops here as we want to loop over z cartoon terms included here
for (int j=0; j<3; j++){
for (int k=0; k<3; k++){
for (int m=0; m<3; m++){
for (int n=0; n<3; n++){

    rhs.Ex += epsLLL[0][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_E_term[n][m];
    rhs.Ey += epsLLL[1][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_E_term[n][m];
    rhs.Ez += epsLLL[2][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_E_term[n][m];
    rhs.Bx += epsLLL[0][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_B_term[n][m];
    rhs.By += epsLLL[1][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_B_term[n][m];
    rhs.Bz += epsLLL[2][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_B_term[n][m];
}}}}

// were done, phew!

}

#endif /* CCZ4CARTOON_IMPL_HPP_ */
