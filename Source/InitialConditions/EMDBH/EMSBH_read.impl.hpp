/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EMSBH_READ_HPP_)
#error "This file should only be included through EMSBH_read.hpp"
#endif

#ifndef EMSBH_READ_IMPL_HPP_
#define EMSBH_READ_IMPL_HPP_

#include "EMSBHSolution_read.hpp" //for EMSBHSolution class

inline EMSBH_read::EMSBH_read(EMSBH_params_t a_params_EMSBH,
                    CouplingFunction::params_t a_params_coupling_function,
                            double a_G_Newton, double a_dx, int a_verbosity)
    : m_dx(a_dx), m_G_Newton(a_G_Newton),
      m_params_EMSBH(a_params_EMSBH),
      m_params_coupling_function(a_params_coupling_function), m_verbosity(a_verbosity)
{
    m_data_path = m_params_EMSBH.data_path;
}

void EMSBH_read::compute_1d_solution(const double max_r)
{
    try
    {
        // set initial parameters and then run the solver
        pout() << "run m_1d_sol.main()" << endl;
        m_1d_sol.main(m_data_path);
        pout() << "completed m_1d_sol.main()" << endl;
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t> void EMSBH_read::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<EinsteinMaxwellDilatonField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_EMSBH.star_centre);


    // binary parameters (can be NO binary too)
    bool binary = m_params_EMSBH.binary;
    double separation = m_params_EMSBH.separation;

    // matter conversion from unit conventions
    double root_kappa = 1./sqrt(8.*M_PI);

    // the kroneka delta
    double kroneka_delta[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};



    ////////////////////////////
    // read 1st black hole
    ////////////////////////////

    // coord objects
    double x = coords.x - 0.5 * separation;
    double z = 0. // coords.z;
    double y = coords.y;
    double cart_coords[3] = {x, y, z};

    // radii and safe (divisible) radii
    double r = sqrt(x * x + y * y + z * z);
    double safe_r = sqrt(x * x + y * y + z * z + 10e-20);
    double rho = sqrt(x * x + y * y);
    double safe_rho = sqrt(x * x + y * y + 10e-20);

    // trig functions
    double sintheta = rho/safe_r;
    double costheta = z/safe_r;
    double sinphi = y/safe_rho;
    double cosphi = x/safe_rho;

    // jacobeans
    double dx_dr = cosphi*sintheta;
    double dy_dr = sinphi*sintheta;
    double dz_dr = costheta;
    double dr_dx = (x / safe_r);
    double dr_dy = (y / safe_r);
    double dr_dz = (z / safe_r);

    // partial cartesian coords (i) by partial polar coords (j)
    // dxc_dxp[i][j]
    // i = {x,y,z}
    // j = {r,th,ph}
    double dxc_dxp[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    dxc_dxp[0][0] = dx_dr;
    dxc_dxp[1][0] = dy_dr;
    dxc_dxp[2][0] = dz_dr;
    dxc_dxp[0][1] = r * cosphi * costheta;
    dxc_dxp[1][1] = r * sinphi * costheta;
    dxc_dxp[2][1] = -r * sintheta;
    dxc_dxp[0][2] = -r * sinphi * sintheta;
    dxc_dxp[1][2] = r * cosphi * sintheta;
    dxc_dxp[2][2] = 0.; // dz/dphi=0

    // partial polar coords (i) by partial cartesian coords (j)
    // dxp_dxc[i][j]
    // i = {r,th,ph}
    // j = {x,y,z}
    double dxp_dxc[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    dxp_dxc[0][0] = dr_dx;
    dxp_dxc[0][1] = dr_dy;
    dxp_dxc[0][2] = dr_dz;
    dxp_dxc[1][0] = x * z / (safe_r * safe_r * safe_rho);
    dxp_dxc[1][1] = y * z / (safe_r * safe_r * safe_rho);
    dxp_dxc[1][2] = - rho / (safe_r * safe_r);
    dxp_dxc[2][0] = - y / (safe_rho * safe_rho);
    dxp_dxc[2][1] = x / (safe_rho * safe_rho);
    dxp_dxc[2][2] = 0.; // dphi/dz=0

    // gamma_polar = 1/X^2 (a dr^2 + b r^2 (dth^2 + sin^2(th) dph^2))
    double X = m_1d_sol.get_value_interp(m_1d_sol.X,r);
    double a = m_1d_sol.get_value_interp(m_1d_sol.a,r);
    double b = m_1d_sol.get_value_interp(m_1d_sol.b,r);

    // spatial metrics
    double gamma_polar[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    gamma_polar[0][0] = a/X/X;
    gamma_polar[1][1] = b*r*r/X/X;
    gamma_polar[2][2] = b*pow(sintheta*r/X,2);
    double gamma_polar_inv[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    gamma_polar_inv[0][0] = 1./gamma_polar[0][0];
    gamma_polar_inv[1][1] = 1./gamma_polar[1][1];
    gamma_polar_inv[2][2] = 1./gamma_polar[2][2];
    double gamma[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gamma_inv[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

    // create cartesian metric from transformation of polar version
    // explicit loops here include z direction
    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
    for (int m=0; m<3; m++){
    for (int n=0; n<3; n++)
    {
        gamma[i][j] += gamma_polar[m][n]*dxp_dxc[m][i]*dxp_dxc[n][j];
        //gamma_inv[i][j] += gamma_polar_inv[m][n]*dxc_dxp[m][i]*dxc_dxp[n][j];
        gamma_inv[i][j] += gamma_polar_inv[m][n]*dxc_dxp[i][m]*dxc_dxp[j][n];
    }}}}}

    // loading the upstairs E^r then lower with gamma_rr = a/(X*X)
    double E_r = m_1d_sol.get_value_interp(m_1d_sol.Er,r)*root_kappa;
    E_r = E_r * gamma_polar[0][0];

    // fabrizio's mixed conformal traceless curvature, then make downstairs
    // THIS ASSUMES DIAGONAL METRIC
    // angular parts must give trace 0
    double AaaUL = m_1d_sol.get_value_interp(m_1d_sol.Aa,r);
    double Aij_polar[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    // this is NOT the conformal A, its simply (K_ij-K gamma_ij/3)
    Aij_polar[0][0] = AaaUL*gamma_polar[0][0];
    Aij_polar[1][1] = -0.5*AaaUL*gamma_polar[1][1];
    Aij_polar[2][2] = -0.5*AaaUL*gamma_polar[2][2];

    // upstairs radial shift
    double betaR = m_1d_sol.get_value_interp(m_1d_sol.shift,r);

    // set gauge vars
    vars.lapse = m_1d_sol.get_value_interp(m_1d_sol.lapse,r);
    vars.shift[0] = dx_dr * betaR;
    vars.shift[1] = dy_dr * betaR;
    // no z shift
    // vars.shift[2] = dz_dr * betaR;

    // set geometry
    // metric det (dummy_chi)
    double dummy_chi = gamma[0][0]*gamma[1][1]*gamma[2][2];
    dummy_chi += 2.*gamma[0][2]*gamma[0][1]*gamma[1][2];
    dummy_chi -= gamma[0][0]*pow(gamma[1][2],2);
    dummy_chi -= gamma[1][1]*pow(gamma[2][0],2);
    dummy_chi -= gamma[2][2]*pow(gamma[0][1],2);
    // assign real chi from power 1/3 of metric det
    vars.chi = pow(dummy_chi,-1./3.);
    vars.K = m_1d_sol.get_value_interp(m_1d_sol.K,r);

    // don't need to loop z direction
    FOR2(i,j)
    {
        vars.h[i][j] = vars.chi*gamma[i][j];
        vars.A[i][j] = 0.;
    }
    //cartoon terms
    vars.hww = vars.chi*gamma[2][2];
    vars.Aww = 0.;

    // explicit loops here include z direction, note i,j, dont include z
    for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
    for (int m=0; m<3; m++){
    for (int n=0; n<3; n++)
    {
        // conformal decomposition here with chi
        vars.A[i][j] += vars.chi * dxp_dxc[m][i] * dxp_dxc[n][j] * Aij_polar[m][n];
    }
    for (int m=0; m<3; m++){
    for (int n=0; n<3; n++)
    {
        // cartoon term
        vars.Aww += vars.chi * dxp_dxc[m][2] * dxp_dxc[n][2] * Aij_polar[m][n];
    }

    // scalar field
    vars.phi = m_1d_sol.get_value_interp(m_1d_sol.phi,r)*root_kappa;
    // minus sign to match Fabrizio's convention
    vars.Pi = -m_1d_sol.get_value_interp(m_1d_sol.pi,r)*root_kappa;

    // electric field
    vars.Ex = dr_dx * E_r;
    vars.Ey = dr_dy * E_r;
    vars.Ez = dr_dz * E_r;


    ////////////////////////////
    // zeros (magnetics field and electromag constrinats)
    ////////////////////////////

    vars.Xi = 0.; // magnetic constraint var
    vars.Lambda = 0.; // leccy constraint var

    // magnetic fields
    vars.Bx = 0.;
    vars.By = 0.;
    vars.Bz = 0.;



    ////////////////////////////
    // flat space spherical harmonic perturbation
    ////////////////////////////

    // // reset coords
    // x = coords.x;
    // z = 0. // coords.z;
    // y = coords.y;
    // r = sqrt(x * x + y * y + z * z);
    // safe_r = sqrt(x * x + y * y + z * z + 10e-20);
    //
    // // new z-plane radius
    // rho = sqrt(x * x + y * y);
    // safe_rho = sqrt(x * x + y * y + 10e-20);
    //
    // // trig
    // sintheta = rho/safe_r;
    // costheta = z/safe_r;
    // sinphi = y/safe_rho;
    // cosphi = x/safe_rho;
    // double cos2phi = cosphi * cosphi - sinphi * sinphi;
    //
    // //spherial harmonic
    // // Y22 is actually Y22 + Y2-2 halved
    // double Y20 = 0.25 * sqrt(5./M_PI) * (3. * costheta * costheta - 1.);
    // double Y22 = 0.25 * sqrt(7.5/M_PI) * (cos2phi * sintheta * sintheta);
    //
    // // inward travelling thin shell params
    // // amplitude is A/r in flat space
    // double A = m_params_EMSBH.Ylm_amplitude;
    // // standard deviation, width of gaussian shell
    // double sig = m_params_EMSBH.Ylm_thickness;
    // // initial radiause of shell
    // double r_0 = m_params_EMSBH.Ylm_r0;
    //
    // double psi = (A/safe_r) * exp(-(r-r_0)*(r-r_0)/(2.*sig*sig));
    //
    // vars.phi += psi;
    // // boosted inwards with Pi - flat space approx no outgoing wave
    // vars.Pi += psi * (r-r_0) / (sig * sig);



    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* EMSBH_READ_IMPL_HPP_ */
