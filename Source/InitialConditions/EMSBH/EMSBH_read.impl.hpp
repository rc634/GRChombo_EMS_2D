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
#include "TensorAlgebra.hpp"

inline EMSBH_read::EMSBH_read(EMSBH_params_t a_params_EMSBH,
                    CouplingFunction::params_t a_params_coupling_function,
                            double a_G_Newton, double a_dx, int a_verbosity)
    : m_dx(a_dx), m_G_Newton(a_G_Newton),
      m_params_EMSBH(a_params_EMSBH),
      m_params_coupling_function(a_params_coupling_function), m_verbosity(a_verbosity)
{
    m_data_path = m_params_EMSBH.data_path;
}

void EMSBH_read::compute_1d_solution()
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

// Compute the value of first bh the initial vars on the grid
template <class data_t> void EMSBH_read::compute(Cell<data_t> current_cell) const
{
    if (m_params_EMSBH.boosted)
    {
      // one boost
      compute_boost(current_cell, -1.);
      // other boost
      if (m_params_EMSBH.binary)
      {
          compute_boost(current_cell, 1.);
      }
    }
    else {
        if (m_params_EMSBH.binary)
        {
          compute1(current_cell); // first bh
          compute2(current_cell); // second bh
          compute3(current_cell); // superposition fixes/fudges
        }
        else
        {
          compute_single(current_cell); // pure solution for one bh
        }
    }
}

// Compute the value of first bh the initial vars on the grid
template <class data_t> void EMSBH_read::compute1(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    // Load variables (should be set to zero)
    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params_EMSBH.star_centre);


    // binary parameters
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
    double z = 0.; // coords.z;
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
    for (int n=0; n<3; n++){
        gamma[i][j] += gamma_polar[m][n]*dxp_dxc[m][i]*dxp_dxc[n][j];
        gamma_inv[i][j] += gamma_polar_inv[m][n]*dxc_dxp[i][m]*dxc_dxp[j][n];
    }}}}

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
    //vars.lapse = m_1d_sol.get_value_interp(m_1d_sol.lapse,r);
    vars.shift[0] = dx_dr * betaR;
    vars.shift[1] = dy_dr * betaR;
    // no z shift

    // curvature
    vars.K = m_1d_sol.get_value_interp(m_1d_sol.K,r);

    // don't loop z direction
    FOR2(i,j)
    {
        vars.h[i][j] = gamma[i][j];
        vars.A[i][j] = 0.;
    }
    //cartoon terms
    vars.hww = gamma[2][2];
    vars.Aww = 0.;

    // explicit loops here include z direction, note i,j, dont include z
    for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
    for (int m=0; m<3; m++){
    for (int n=0; n<3; n++){
        // conformal decomposition here with chi done in copute 3
        vars.A[i][j] += dxp_dxc[m][i] * dxp_dxc[n][j] * Aij_polar[m][n];
    }}}}
    for (int m=0; m<3; m++){
    for (int n=0; n<3; n++){
        // cartoon term
        vars.Aww += dxp_dxc[m][2] * dxp_dxc[n][2] * Aij_polar[m][n];
    }}

    // scalar field
    vars.phi = m_1d_sol.get_value_interp(m_1d_sol.phi,r)*root_kappa;
    // minus sign to match Fabrizio's convention
    vars.Pi = -m_1d_sol.get_value_interp(m_1d_sol.pi,r)*root_kappa;

    // electric field
    vars.Ex = dr_dx * E_r;
    vars.Ey = dr_dy * E_r;
    vars.Ez = dr_dz * E_r;

    current_cell.store_vars(vars);
}

// Compute the value of the initial vars on the grid
template <class data_t> void EMSBH_read::compute2(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    // Load variables (should be set to zero)
    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params_EMSBH.star_centre);


    // binary parameters (can be NO binary too)
    bool binary = m_params_EMSBH.binary;
    if (binary) {
    double separation = m_params_EMSBH.separation;

    // matter conversion from unit conventions
    double root_kappa = 1./sqrt(8.*M_PI);

    // the kroneka delta
    double kroneka_delta[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};



    ////////////////////////////
    // read 2st black hole (if it exists), no boosts
    ////////////////////////////

    // coord objects
    double x = coords.x + 0.5 * separation;
    double z = 0.; // coords.z;
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
    for (int n=0; n<3; n++){
        gamma[i][j] += gamma_polar[m][n]*dxp_dxc[m][i]*dxp_dxc[n][j];
        gamma_inv[i][j] += gamma_polar_inv[m][n]*dxc_dxp[i][m]*dxc_dxp[j][n];
    }}}}

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
    vars.shift[0] += dx_dr * betaR;
    vars.shift[1] += dy_dr * betaR;
    // no z shift

    // curvature
    vars.K += m_1d_sol.get_value_interp(m_1d_sol.K,r);


    // don't loop z direction
    FOR2(i,j)
    {
        vars.h[i][j] += gamma[i][j] - kroneka_delta[i][j];// remove 1 for superposition only
    }
    //cartoon terms
    vars.hww += gamma[2][2]-1.;;


    // explicit loops here include z direction, note i,j, dont include z
    for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
    for (int m=0; m<3; m++){
    for (int n=0; n<3; n++){
        // conformal decomposition here with chi done in compute 3
        vars.A[i][j] += dxp_dxc[m][i] * dxp_dxc[n][j] * Aij_polar[m][n];
    }}}}
    for (int m=0; m<3; m++){
    for (int n=0; n<3; n++){
        // cartoon term
        vars.Aww += dxp_dxc[m][2] * dxp_dxc[n][2] * Aij_polar[m][n];
    }}

    // scalar field
    vars.phi += m_1d_sol.get_value_interp(m_1d_sol.phi,r)*root_kappa;
    // minus sign to match Fabrizio's convention
    vars.Pi += -m_1d_sol.get_value_interp(m_1d_sol.pi,r)*root_kappa;

    // electric field
    vars.Ex += dr_dx * E_r;
    vars.Ey += dr_dy * E_r;
    vars.Ez += dr_dz * E_r;


    ////////////////////////////
    // zeros (magnetics field and electromag constrinats)
    ////////////////////////////

    // Xi, Lambda, B_i alrady 0

    current_cell.store_vars(vars);

}}

// Compute the value of the initial vars on the grid
template <class data_t> void EMSBH_read::compute3(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    // Load variables (should be set to zero)
    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params_EMSBH.star_centre);


    ////////////////////////////
    // generic stuff
    ////////////////////////////

    vars.Xi = 0.; // magnetic constraint var
    vars.Lambda = 0.; // leccy constraint var

    // magnetic fields
    vars.Bx = 0.;
    vars.By = 0.;
    vars.Bz = 0.;


    /////////////////
    // second adjustment of conformal factor from superposition errors
    /////////////////
    double det_h = (vars.h[0][0]*vars.h[1][1]-vars.h[1][0]*vars.h[0][1])*vars.hww;

    vars.chi = pow(det_h,-1./3.);
    vars.lapse = sqrt(vars.chi);


    ////////////////////////////
    // fix bullshit from superposition
    ////////////////////////////

    FOR2(i,j)
    {
        vars.h[i][j] *= vars.chi;
        vars.A[i][j] *= vars.chi;
    }
    vars.Aww *= vars.chi;
    vars.hww *= vars.chi;

    vars.lapse = sqrt(vars.chi);

    current_cell.store_vars(vars);

}







//////////////

/////////////

///// working old version down here

/////////////

////////////






// Compute the value of first bh the initial vars on the grid
template <class data_t> void EMSBH_read::compute_single(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    // Load variables (should be set to zero)
    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params_EMSBH.star_centre);


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
    double z = 0.; // coords.z;
    double y = coords.y;
    double cart_coords[3] = {x, y, z};

    // radii and safe (divisible) radii
    double epsilon_small = 10e-20;
    double r = sqrt(x * x + y * y + z * z);
    double safe_r = sqrt(x * x + y * y + z * z + epsilon_small);
    double rho = sqrt(x * x + y * y);
    double safe_rho = sqrt(x * x + y * y + epsilon_small);

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
    for (int n=0; n<3; n++){
        gamma[i][j] += gamma_polar[m][n]*dxp_dxc[m][i]*dxp_dxc[n][j];
        gamma_inv[i][j] += gamma_polar_inv[m][n]*dxc_dxp[i][m]*dxc_dxp[j][n];
    }}}}

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
    // metric det = dummy_chi
    double dummy_chi = gamma[0][0]*gamma[1][1]*gamma[2][2];
    dummy_chi += 2.*gamma[0][2]*gamma[0][1]*gamma[1][2];
    dummy_chi -= gamma[0][0]*pow(gamma[1][2],2);
    dummy_chi -= gamma[1][1]*pow(gamma[2][0],2);
    dummy_chi -= gamma[2][2]*pow(gamma[0][1],2);
    // assign real chi from power -1/3 of metric det
    vars.chi = pow(dummy_chi,-1./3.);
    vars.K = m_1d_sol.get_value_interp(m_1d_sol.K,r);

    // manual gauge forcing - probably optional for one object
    vars.lapse = sqrt(vars.chi);
    // vars.shift[0] = 0.;
    // vars.shift[1] = 0.;

    // don't loop z direction
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
    for (int n=0; n<3; n++){
        // conformal decomposition here with chi done in copute 3
        vars.A[i][j] += vars.chi * dxp_dxc[m][i] * dxp_dxc[n][j] * Aij_polar[m][n];
    }}}}
    for (int m=0; m<3; m++){
    for (int n=0; n<3; n++){
        // cartoon term
        vars.Aww += vars.chi * dxp_dxc[m][2] * dxp_dxc[n][2] * Aij_polar[m][n];
    }}

    // scalar field
    vars.phi = m_1d_sol.get_value_interp(m_1d_sol.phi,r)*root_kappa;
    // minus sign to match Fabrizio's convention
    vars.Pi = -m_1d_sol.get_value_interp(m_1d_sol.pi,r)*root_kappa;

    // electric field
    vars.Ex = dr_dx * E_r;
    vars.Ey = dr_dy * E_r;
    vars.Ez = dr_dz * E_r;



    current_cell.store_vars(vars);
}



///////////////////////////////////////

//     BOOSTED

//     INITIAL

//      DATA

///////////////////////////////////////



// Compute the value of first bh the initial vars on the grid
template <class data_t> void EMSBH_read::compute_boost(Cell<data_t> current_cell
                                                          , double a_sign) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    // Load variables (should be set to zero)
    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params_EMSBH.star_centre);


    // binary parameters
    bool binary = m_params_EMSBH.binary;
    double separation = m_params_EMSBH.separation;
    double rapidity = m_params_EMSBH.rapidity;

    // matter conversion from unit conventions
    const double root_kappa = 1./sqrt(8.*M_PI);

    // the kroneka delta
    const double kroneka_delta[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};


    // protected single letter names
    // a, b, r, x, y, z, X



    ////////////////////////////
    // 4 - METRIC
    ////////////////////////////

    // 4-metrics for boosting numerically
    double gmunu[4][4] = {
                 {0., 0., 0., 0.},
                 {0., 0., 0., 0.},
                 {0., 0., 0., 0.},
                 {0., 0., 0., 0.}
    };
    double gmunu_boosted[4][4] = {
                {0., 0., 0., 0.},
                {0., 0., 0., 0.},
                {0., 0., 0., 0.},
                {0., 0., 0., 0.}
    };
    double gmunu_inv[4][4] = {
                 {0., 0., 0., 0.},
                 {0., 0., 0., 0.},
                 {0., 0., 0., 0.},
                 {0., 0., 0., 0.}
    };
    double gmunu_inv_boosted[4][4] = {
                {0., 0., 0., 0.},
                {0., 0., 0., 0.},
                {0., 0., 0., 0.},
                {0., 0., 0., 0.}
    };

    // remember ADM decomp
    // g_tt = -lapse^2 + shift^2
    // g_ti = shift_i
    // g_ij = gamma_ij




    ////////////////////////////
    // 4 - BOOST MATRICES
    ////////////////////////////

    // boost matrix
    const double c_ = cosh(rapidity);
    // use sign = 1 to denote positive velocity boost
    // aka BH moves towards the right
    const double s_ = sinh(a_sign * rapidity);

    // Boost matrix
    // remember time index is 3.
    // For co-tensor components
    // from rest frame to lab frame
    const double Lambda[4][4] = {
                {c_, 0., 0., -s_},
                {0., 1., 0., 0.},
                {0., 0., 1., 0.},
                {-s_, 0., 0., c_}
    };
    // inverse boosts upsairs indices
    const double Lambda_inv[4][4] = {
                {c_, 0., 0., s_},
                {0., 1., 0., 0.},
                {0., 0., 1., 0.},
                {s_, 0., 0., c_}
    };

    // identity matrix for shits and giggles
    const double I4[4][4] = {
                {1., 0., 0., 0.},
                {0., 1., 0., 0.},
                {0., 0., 1., 0.},
                {0., 0., 0., 1.}
    };



    ////////////////////////////
    // COORDINATES
    ////////////////////////////


    // coord objects
    // sign = 1 corresponds to bh positioned to the left of the binary
    double x = c_ * (coords.x + a_sign * 0.5 * separation);
    double z = 0.; // coords.z;
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

    double cosphi2 = cosphi*cosphi;
    double costheta2 = costheta*costheta;
    double sinphi2 = sinphi*sinphi;
    double sintheta2 = sintheta*sintheta;

    // jacobeans
    double dx_dr = cosphi*sintheta;
    double dy_dr = sinphi*sintheta;
    double dz_dr = costheta;
    double dr_dx = (x / safe_r);
    double dr_dy = (y / safe_r);
    double dr_dz = (z / safe_r);
    double dth_dx = x*z/(safe_r*safe_r*safe_rho);
    double dth_dy = y*z/(safe_r*safe_r*safe_rho);
    double dth_dz = -rho/(safe_r*safe_r);
    double dph_dx = -y/(safe_rho*safe_rho);
    double dph_dy = x/(safe_rho*safe_rho);
    double dph_dz = 0.;

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



    ////////////////////////////
    // LOAD RADIAL SOLUTION
    ////////////////////////////

    // gamma_polar = 1/X^2 (a dr^2 + b r^2 (dth^2 + sin^2(th) dph^2))
    double X = m_1d_sol.get_value_interp_o4(m_1d_sol.X,r);
    double a = m_1d_sol.get_value_interp_o4(m_1d_sol.a,r);
    double b = m_1d_sol.get_value_interp_o4(m_1d_sol.b,r);
    double alpha = m_1d_sol.get_value_interp_o4(m_1d_sol.lapse,r);
    double betaR = m_1d_sol.get_value_interp_o4(m_1d_sol.shift,r);
    double dXdr = m_1d_sol.get_deriv_interp_o4(m_1d_sol.X,r);
    double dadr = m_1d_sol.get_deriv_interp_o4(m_1d_sol.a,r);
    double dbdr = m_1d_sol.get_deriv_interp_o4(m_1d_sol.b,r);
    double dalpha_dr = m_1d_sol.get_deriv_interp_o4(m_1d_sol.lapse,r);
    double dbetaR_dr = m_1d_sol.get_deriv_interp_o4(m_1d_sol.shift,r);
    //
    double X2 = X*X;
    double X3 = X2*X;
    double s2 = s_*s_;
    double c2 = c_*c_;
    double a2 = alpha*alpha;
    // some time derivs - from gauge conditions of Fabrizio
    double dalpha_dt = 0.;
    double dbetaR_dt = m_1d_sol.get_value_interp_o4(m_1d_sol.Br,r);
    double d_gamma_rr_dr = dadr/X2 - 2.*a*dXdr/X3;
    double dt_shift[3] = {x/safe_r*dbetaR_dt,
                          y/safe_r*dbetaR_dt,
                          z/safe_r*dbetaR_dt};



    //////////////////////////////////
    // CONSTRUCT SPATIAL METRICS
    //////////////////////////////////

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
    for (int n=0; n<3; n++){
        gamma[i][j] += gamma_polar[m][n]*dxp_dxc[m][i]*dxp_dxc[n][j];
        gamma_inv[i][j] += gamma_polar_inv[m][n]*dxc_dxp[i][m]*dxc_dxp[j][n];
    }}}}
    // spatial volume element for B boost field
    double root_gamma = 0.;
    // this is the same
    // root_gamma += gamma[0][0] * (gamma[1][1]*gamma[2][2]
    //                           - gamma[1][2]*gamma[2][1]);
    // root_gamma += gamma[0][1] * (gamma[1][2]*gamma[2][0]
    //                           - gamma[1][0]*gamma[2][2]);
    // root_gamma += gamma[0][2] * (gamma[1][0]*gamma[2][1]
    //                           -gamma[1][1]*gamma[2][0]);
    root_gamma = sqrt(a)*b/X3;



    //////////////////////////////////////////////
    // Lie derivative of gamma_ij w.r.t. shift
    //////////////////////////////////////////////
    double LieBeta_gamma_polar[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double LieBeta_gamma_cart[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

    // total lie derivative
    LieBeta_gamma_polar[0][0] = (X*betaR*dadr - 2.*a*betaR*dXdr + 2.*a*X*dbetaR_dr)/X3;
    LieBeta_gamma_polar[1][1] = r*betaR*(r*X*dbdr + 2.*b*(X-r*dXdr))/X3;
    LieBeta_gamma_polar[2][2] = r*sintheta2*betaR*(r*X*dbdr + 2.*b*(X-r*dXdr))/X3;

    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
    for (int m=0; m<3; m++){
    for (int n=0; n<3; n++){
       LieBeta_gamma_cart[i][j] += dxp_dxc[n][j] * dxp_dxc[m][i] * LieBeta_gamma_polar[m][n];
    }}}}



    ////////////////////////////////
    // POPULATE 4 METRIC
    ////////////////////////////////

    // note that time index is 3 not 0!
    // x,y,z,t = 0,1,2,3 -- throughout !

    // cartesian shift setup

    double beta_squared = 0.;
    double beta_U[3] = {0.,0.,0.}; // upper components
    double beta_L[3] = {0.,0.,0.}; // lower components
    beta_U[0] = dx_dr * betaR;
    beta_U[1] = dy_dr * betaR;
    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
        beta_L[i] += gamma[i][j] * beta_U[j];
    }}
    for (int i=0; i<3; i++){
        beta_squared += beta_L[i] * beta_U[i];
    }

    // ADM metric formula for creating 4d metric from 3d component
    // CAREFUL, TIME INDEX IS 3 HERE !!!!!!!

    // g_tt, time-time part
    gmunu[3][3] = - alpha * alpha + beta_squared;
    gmunu_inv[3][3] = -1./(alpha*alpha);

    // g_it, space-time mixed part
    for (int i=0; i<3; i++){
        gmunu[3][i] = beta_L[i];
        gmunu[i][3] = beta_L[i];
        gmunu_inv[3][i] = beta_U[i] / (alpha*alpha);
        gmunu_inv[i][3] = beta_U[i] / (alpha*alpha);
    }

    // g_ij, space space part
    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
        gmunu[i][j] = gamma[i][j];
        gmunu_inv[i][j] = gamma_inv[i][j] - beta_U[i]*beta_U[j]/(alpha*alpha);
    }}



    ///////////////////////////////
    // DO THE BOOST ON 4-METRIC
    ///////////////////////////////

    // Do the boost ! // need 0-3 indices here, hence <4
    for (int m=0; m<4; m++){
    for (int n=0; n<4; n++){
    for (int t=0; t<4; t++){
    for (int s=0; s<4; s++){
        gmunu_boosted[m][n] += Lambda[m][t] * Lambda[n][s]
                             * gmunu[t][s];
        gmunu_inv_boosted[m][n] += Lambda_inv[m][t] * Lambda_inv[n][s]
                                 * gmunu_inv[t][s];
    }}}}

    // no boost 1
    // for (int m=0; m<4; m++){
    // for (int n=0; n<4; n++){
    // for (int q=0; q<4; q++){
    // for (int s=0; s<4; s++){
    //     gmunu_boosted[m][n] += I4[m][q] * I4[n][s] * gmunu[q][s];
    //     gmunu_inv_boosted[m][n] += I4[m][q] * I4[n][s] * gmunu_inv[q][s];
    // }}}}

    // no boost 2
    // for (int m=0; m<4; m++){
    // for (int n=0; n<4; n++){
    //     gmunu_boosted[m][n] = gmunu[m][n];
    //     gmunu_inv_boosted[m][n] = gmunu_inv[m][n];
    // }}



    ///////////////////////////////////////
    // ADM DECOMPOSE BOOSTED 4-METRIC
    ///////////////////////////////////////

    // 3+1 ADM decomp to get boosted lapse, shift, gamma_ij

    // the boosted shift
    double boosted_beta_U[3] = {0.,0.,0.}; // upper components
    double boosted_beta_L[3] = {0.,0.,0.}; // lower components
    Tensor<2, data_t, 3> boosted_gamma_LL = {0.};
    double beta_boost_sqr = 0.;

    for (int i=0; i<3; i++){
        boosted_beta_L[i] = gmunu_boosted[3][i];
    }

    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
        boosted_gamma_LL[i][j] = gmunu_boosted[i][j];
    }}
    // boosted spatial volume element for B boost field
    double root_gamma_tilde = 0.;
    root_gamma_tilde += boosted_gamma_LL[0][0] * (
                              boosted_gamma_LL[1][1]*boosted_gamma_LL[2][2]
                              - boosted_gamma_LL[1][2]*boosted_gamma_LL[2][1]
                                                 );
    root_gamma_tilde += boosted_gamma_LL[0][1] * (
                              boosted_gamma_LL[1][2]*boosted_gamma_LL[2][0]
                              - boosted_gamma_LL[1][0]*boosted_gamma_LL[2][2]
                                                );
    root_gamma_tilde += boosted_gamma_LL[0][2] * (
                              boosted_gamma_LL[1][0]*boosted_gamma_LL[2][1]
                              - boosted_gamma_LL[1][1]*boosted_gamma_LL[2][0]
                                                 );
    root_gamma_tilde = sqrt(root_gamma_tilde);

    // auto boosted_gamma_UU = TensorAlgebra::compute_inverse_sym(boosted_gamma_LL);

    // for (int i=0; i<3; i++){
    // for (int j=0; j<3; j++){
    //     boosted_beta_U[i] += boosted_beta_L[j] * boosted_gamma_UU[j][i];
    // }}
    for (int i=0; i<3; i++){
        boosted_beta_U[i] = -gmunu_inv_boosted[3][i]/gmunu_inv_boosted[3][3];
    }

    for (int i=0; i<3; i++){
        beta_boost_sqr += boosted_beta_L[i] * boosted_beta_U[i];
    }

    // needed to calculate normal vector and Pi, Kij, Fmunu ...
    double boost_lapse = 0., blt = 0., lapse_diff_err = 0., shiftsqrdiff = 0.;
    double metric_diff = 0.;
    boost_lapse = sqrt(beta_boost_sqr - gmunu_boosted[3][3]);
    blt = sqrt(-1./gmunu_inv_boosted[3][3]); // boost lapse test
    lapse_diff_err = abs(blt-alpha);
    shiftsqrdiff = abs(beta_boost_sqr-beta_squared);

    for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
        metric_diff += pow(gmunu_boosted[i][j]-gmunu[i][j],2);
    }}




    ///////////////////////////////////////
    // STORE BOOSTED SHIFT AND GAMMA_IJ
    ///////////////////////////////////////

    // lapse likely not used and later over-written
    vars.lapse += alpha;
    vars.shift[0] += boosted_beta_U[0];
    vars.shift[1] += boosted_beta_U[1];
    // no z shift

    // test evaluation of shift y
    // double test_beta_y = boost_lapse*boost_lapse*(c_/alpha/alpha*beta_U[1]);
    // test_beta_y += boost_lapse*boost_lapse*(s_*gmunu_inv[0][1]);
    // vars.shift[1] += -test_beta_y;

    // don't loop z direction
    // store physical metric for now
    FOR2(i,j)
    {
        vars.h[i][j] += gmunu_boosted[i][j]; // physical metric for now
    }
    //cartoon terms
    vars.hww += gmunu_boosted[2][2]; // physical metric for now




    ///////////////////////////
    // SCALAR FIELD
    //////////////////////////

    vars.phi += m_1d_sol.get_value_interp_o4(m_1d_sol.phi,r)*root_kappa;



    // derivatives are rest frame
    double dphi_dr = m_1d_sol.get_deriv_interp_o4(m_1d_sol.phi,r)*root_kappa;
    // minus sign from different convention to Fabrizio
    double default_pi = -m_1d_sol.get_value_interp_o4(m_1d_sol.pi,r)*root_kappa;
    double dphi_dx = dphi_dr * dr_dx;
    double dphi_dy = dphi_dr * dr_dy;

    // close to zero but might not quite be
    double dphi_dt = 0.;
    dphi_dt = - alpha * default_pi + beta_U[0]*dphi_dx + beta_U[1]*dphi_dy;

    // boost done implicitly here
    vars.Pi += (1./boost_lapse) * (
                    - (c_ + boosted_beta_U[0] * s_) * dphi_dt
                    + (s_ + boosted_beta_U[0] * c_) * dphi_dx
                          + boosted_beta_U[1]       * dphi_dy);

    // hard remove singularity in Pi
    // if (vars.Pi > 0.018)
    // {
    //     vars.Pi = -0.019;
    // }

    // used this to check if the calculated Pi is similar to the read one
    // it is almost identical
    // only uncomment for testing
    // vars.Pi -= default_pi;




    ////////////////////////
    // ELECTROMAGNETISM
    ////////////////////////

    // loading the upstairs E^r then lower with gamma_rr = a/(X*X)
    double EUr = m_1d_sol.get_value_interp_o4(m_1d_sol.Er,r)*root_kappa;
    double E_r = EUr * gamma_polar[0][0];

    // only included non-zero temrs in boost, this is not generic!!
    double Ex = dr_dx * E_r;
    double Ey = dr_dy * E_r;
    double Ez = dr_dz * E_r; // shoudl be 0

    // boost calcualted by hand (differential forms make this easier)
    // vars.Ex += Ex * alpha / boost_lapse;
    // vars.Ey += c_ * Ey * alpha / boost_lapse;
    // vars.Ez += c_ * Ez * alpha / boost_lapse;
    // vars.Bx += 0.;
    // vars.By += 0.;
    // //vars.Bz += alpha * s_ * Ey * gmunu_boosted[2][2] / root_gamma_tilde;
    // vars.Bz += boost_lapse * s_ * Ey * gmunu_boosted[2][2] / root_gamma;

    vars.Ex += a * cosphi * alpha * EUr / X2 / boost_lapse;
    vars.Ey += a * c_ * sinphi * alpha * EUr / X2 / boost_lapse;
    vars.Bz += a * sinphi * s_ * alpha * EUr * boosted_gamma_LL[2][2] / X2 / root_gamma_tilde;

    // shubbezeros
    vars.Ez += 0.;
    vars.Bx += 0.;


    double shubbezero = 0., dummyvar=0.;

    dummyvar = cosphi2*c_*s_*X2/a + c_*sinphi2*s_*X2/b + (-c_*s_ + cosphi*betaR*(c2+s2-cosphi*c_*s_*betaR))/a2;

    shubbezero = pow(dummyvar*boost_lapse*boost_lapse-boosted_beta_U[0],2);

    vars.By += shubbezero;

    // list of ok'ed components
    // gmunu boost : passed
    // gmunu * gmunuinv = 4 : passed
    // gmunu inverse : seems fine 
    // lapse : passed
    // shift : passed

    // // numerical F boost
    // double faraday[4][4] = {
    //             {0., 0., 0., 0.},
    //             {0., 0., 0., 0.},
    //             {0., 0., 0., 0.},
    //             {0., 0., 0., 0.}
    // };
    // faraday[3][0] = -alpha*Ex;
    // faraday[0][3] = alpha*Ex;
    // faraday[3][1] = -alpha*Ey;
    // faraday[1][3] = alpha*Ey;
    //
    // double faraday_boosted[4][4] = {
    //             {0., 0., 0., 0.},
    //             {0., 0., 0., 0.},
    //             {0., 0., 0., 0.},
    //             {0., 0., 0., 0.}
    // };
    //
    // for (int m=0; m<4; m++){
    // for (int n=0; n<4; n++){
    // for (int r=0; r<4; r++){
    // for (int s=0; s<4; s++){
    //   faraday_boosted[m][n] +=        Lambda[m][r]
    //                                 * Lambda[n][s]
    //                                 * faraday[r][s];
    // }}}}
    //vars.Bz += boosted_gamma_LL[2][2] * faraday_boosted[0][1] / root_gamma_tilde;






    ////////////////////////////////////
    // LOAD EXTRINSIC CURVATURE K_IJ
    ////////////////////////////////////

    // fabrizio's mixed conformal traceless curvature, then make downstairs
    // THIS ASSUMES DIAGONAL METRIC
    // angular parts must give trace 0
    double AaaUL = m_1d_sol.get_value_interp(m_1d_sol.Aa,r);
    double inputK = m_1d_sol.get_value_interp(m_1d_sol.K,r);
    double Kij_polar[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    // this is NOT the conformal A, its simply (K_ij-K gamma_ij/3)
    // becuase i believe that Fabrizios trK is automatically zero
    Kij_polar[0][0] = (AaaUL - inputK/3.)*gamma_polar[0][0];
    Kij_polar[1][1] = (-0.5*AaaUL- inputK/3.)*gamma_polar[1][1];
    Kij_polar[2][2] = (-0.5*AaaUL- inputK/3.)*gamma_polar[2][2];

    // could be useful for debugging
    double Kij_rest_frame[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
    for (int m=0; m<3; m++){
    for (int n=0; n<3; n++){
       Kij_rest_frame[i][j] += dxp_dxc[n][j] * dxp_dxc[m][i] * Kij_polar[m][n];
    }}}}

    double L_beta_gamma_rr = 2. * gamma_polar[0][0] * dbetaR_dr;
    L_beta_gamma_rr += betaR * d_gamma_rr_dr;
    double d_gamma_rr_dt = L_beta_gamma_rr-2.*alpha*Kij_polar[0][0];




    ////////////////////////////////////////
    // BOOSTED EXTRINSIC CURVATURE SETUP
    ////////////////////////////////////////

    // time derivative of gamma_ij and stuff
    double dt_gamma[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double dt_tilde_gamma_boost[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

    // Lie derivative of boosted gamma_ij
    double L_b_gamma_boost[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

    // metric gradients -- big objects ....
    double d_4metric[4][4][4] = {
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}},
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}},
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}},
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}
    };
    double d_4metric_boosted[4][4][4] = {
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}},
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}},
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}},
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}
    };
    double d_inv_4metric_boosted[4][4][4] = {
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}},
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}},
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}},
      {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}
    };
    double dr_4metric[4][4] = {{0., 0., 0., 0.},
                               {0., 0., 0., 0.},
                               {0., 0., 0., 0.},
                               {0., 0., 0., 0.}};
    double dth_4metric[4][4] = {{0., 0., 0., 0.},
                               {0., 0., 0., 0.},
                               {0., 0., 0., 0.},
                               {0., 0., 0., 0.}};
    double dph_4metric[4][4] = {{0., 0., 0., 0.},
                               {0., 0., 0., 0.},
                               {0., 0., 0., 0.},
                               {0., 0., 0., 0.}};
    double dt_4metric[4][4] = {{0., 0., 0., 0.},
                               {0., 0., 0., 0.},
                               {0., 0., 0., 0.},
                               {0., 0., 0., 0.}};


    // time derivs of 4d metric
    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
        dt_4metric[i][j] = LieBeta_gamma_cart[i][j] - 2.*alpha*Kij_rest_frame[i][j];
    }}
    // tt component
    dt_4metric[3][3] = -2.*alpha*dalpha_dt
                             + betaR*betaR*d_gamma_rr_dt
                             + 2.*gamma_polar[0][0]*betaR*dbetaR_dr;
    // CRUCIAL TO ENSURE THE SPATIAL PART OF dt_4metric is already calculated
    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
        dt_4metric[3][j] += gamma_polar[i][j] * dt_shift[i] + beta_U[i]*dt_4metric[i][j];
    }
    }
    // symmetriser
    dt_4metric[0][3] = dt_4metric[3][0];
    dt_4metric[1][3] = dt_4metric[3][1];
    dt_4metric[2][3] = dt_4metric[3][2];



    ////////////////////////////////////////
    // CALCULATE REST FRAME METRIC DERIVS
    ////////////////////////////////////////

    // oh boy here we go .... time to calcualte d_g ...
    // in cartesian gauge, but expressed in polars for ease

    /////////////////////
    // // d g \ d r
    // xx
    dr_4metric[0][0] = (cosphi2*dadr+sinphi2*dbdr)/X2;
    dr_4metric[0][0] += -2.*(a*cosphi2 + b*sinphi2)*dXdr/X3;
    // xy
    dr_4metric[0][1] = cosphi*sinphi*(dadr-dbdr)/X2;
    dr_4metric[0][1] += -2.*cosphi*sinphi*(a-b)*dXdr/X3;
    dr_4metric[1][0] = dr_4metric[0][1];
    // xz = 0
    // yy
    dr_4metric[1][1] = (sinphi2*dadr + cosphi2*dbdr)/X2;
    dr_4metric[1][1] += -2.*(a*sinphi2 + b*cosphi2)*dXdr/X3;
    // yz = 0
    // zz
    dr_4metric[2][2] = dbdr/X2 - 2.*b*dXdr/X3;
    // tt
    dr_4metric[3][3] = -2.*a*betaR*betaR*dXdr/X3 - 2.*alpha*dalpha_dr;
    dr_4metric[3][3] += betaR * (betaR * dadr + 2. * a * dbetaR_dr)/X2;
    // tx
    dr_4metric[0][3] = cosphi*(-2.*a*betaR*dXdr + X*(betaR*dadr + a*dbetaR_dr))/X3;
    dr_4metric[3][0] = dr_4metric[0][3];
    // ty
    dr_4metric[1][3] = sinphi*(-2.*a*betaR*dXdr + X*(betaR*dadr + a*dbetaR_dr))/X3;
    dr_4metric[3][1] = dr_4metric[1][3];
    // tz = 0

    //////////////////////
    // // d g \ d th
    // xx = 0
    // xy = 0
    // xz
    dth_4metric[0][2] = -(a-b)*cosphi/X2;
    dth_4metric[2][0] = dth_4metric[0][2];
    // yy = 0
    // yz
    dth_4metric[1][2] = -(a-b)*sinphi/X2;
    dth_4metric[2][1] = dth_4metric[1][2];
    // zz = 0
    // tt = tx = ty = 0
    // tz
    dth_4metric[2][3] = - a * betaR / X2;
    dth_4metric[3][2] = dth_4metric[2][3];

    //////////////////////
    // // d g \ d ph
    // xx
    dph_4metric[0][0] = -2.*(a-b)*sinphi*cosphi/X2;
    // xy
    dph_4metric[0][1] = (a-b)*(cosphi2-sinphi2)/X2;
    dph_4metric[1][0] = dph_4metric[0][1];
    // xz = 0
    // yy
    dph_4metric[1][1] = 2.*(a-b)*sinphi*cosphi/X2;
    // yz = 0
    // zz = 0
    // tt = tz = 0
    // tx
    dph_4metric[0][3] = - a * sinphi * betaR / X2;
    dph_4metric[3][0] = dph_4metric[0][3];
    // ty
    dph_4metric[1][3] = a * cosphi * betaR / X2;
    dph_4metric[3][1] = dph_4metric[1][3];


    // finally create cartesian 4-metric derivs here
    // loops to 4 here for time!! i,j = 0,1,2,3 <-> x,y,z,t
    for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
        d_4metric[0][i][j] = dr_dx * dr_4metric[i][j]
                           + dth_dx * dth_4metric[i][j]
                           + dph_dx * dph_4metric[i][j];
        d_4metric[1][i][j] = dr_dy * dr_4metric[i][j]
                           + dth_dy * dth_4metric[i][j]
                           + dph_dy * dph_4metric[i][j];
        d_4metric[2][i][j] = dr_dz * dr_4metric[i][j]
                           + dth_dz * dth_4metric[i][j]
                           + dph_dz * dph_4metric[i][j];
        d_4metric[3][i][j] = dt_4metric[i][j];
    }}
    // now boost the results, we have deriv of boosted 4-metric
    // with respect to boosted coords
    for (int m=0; m<4; m++){
    for (int n=0; n<4; n++){
    for (int p=0; p<4; p++){
    for (int q=0; q<4; q++){
    for (int s=0; s<4; s++){
    for (int t=0; t<4; t++){
        d_4metric_boosted[p][m][n] += Lambda[p][t]
                                    * Lambda[m][q]
                                    * Lambda[n][s]
                                    * d_4metric[t][q][s];
    }}}}}}

    ///////// d_inverse_g all LAB frame
    // meaning metric is in LAB gauge and so are derivs
    for (int m=0; m<4; m++){
    for (int n=0; n<4; n++){
    for (int q=0; q<4; q++){
    for (int s=0; s<4; s++){
    for (int t=0; t<4; t++){
        d_inv_4metric_boosted[t][q][s] -= gmunu_inv_boosted[m][q]
                                        * gmunu_inv_boosted[n][s]
                                        * d_4metric_boosted[t][m][n];
    }}}}}

    // calculate the boosted frame shift derivs (quite non trivial)
    double dibj_boosted[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    // d_x beta^x  (all boosted)
    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
        dibj_boosted[i][j] = boost_lapse * boost_lapse *(
                             d_inv_4metric_boosted[i][3][j] // d_i g^{tj}
                             + boosted_beta_U[j] *
                             d_inv_4metric_boosted[i][3][3] // d_i g^{tt}
                                                      );
    }}

    // feed this into dt_tilde_gamma_boost, both gamma and t are LAB frame
    // also calculate the tildebeta dot tildepartial tildegamma of Lie deriv
    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      dt_tilde_gamma_boost[i][j] = d_4metric_boosted[3][i][j];
      L_b_gamma_boost[i][j] =  boosted_beta_U[0] * d_4metric_boosted[0][i][j];
      L_b_gamma_boost[i][j] += boosted_beta_U[1] * d_4metric_boosted[1][i][j];
    }}

    // finish Lie deriv
    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
    for (int k=0; k<3; k++){
      L_b_gamma_boost[i][j] += gmunu_boosted[k][j] * dibj_boosted[i][k];
      L_b_gamma_boost[i][j] += gmunu_boosted[i][k] * dibj_boosted[j][k];
    }}}

    // Finally Kij in LAB frame
    double tilde_Kij[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      tilde_Kij[i][j] = -(dt_tilde_gamma_boost[i][j] - L_b_gamma_boost[i][j])
                         /(2.*boost_lapse);
    }}

    // temporary storage of K_ij in A_ij
    FOR2(i,j)
    {
        vars.A[i][j] += tilde_Kij[i][j];
    }
    //cartoon terms
    vars.Aww += tilde_Kij[2][2];


    //////////////////////////////////////////
    // old stuff left incase needed again
    //////////////////////////////////////////


    // /////////////////////////////////////
    // // shift derivs - rest frame
    double di_betaj[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    //
    // // one way of calculating beta deriv
    // // for (int j=0; j<3; j++){
    // //   di_betaj[0][j] += dr_dx * (dbetaR_dr / safe_r - betaR / pow(safe_r,2)) * cart_coords[j];
    // //   di_betaj[1][j] += dr_dy * (dbetaR_dr / safe_r - betaR / pow(safe_r,2)) * cart_coords[j];
    // //   di_betaj[2][j] += dr_dz * (dbetaR_dr / safe_r - betaR / pow(safe_r,2)) * cart_coords[j];
    // //   // diagonal term
    // //   di_betaj[j][j] += betaR / safe_r;
    // // }
    //
    // // alternative way of calculating the beta deriv
    for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      di_betaj[i][j] += betaR * kroneka_delta[i][j] / safe_r;
      di_betaj[i][j] += cart_coords[i] * cart_coords[j]  *
                              (dbetaR_dr - betaR/safe_r) *
                               pow(safe_r,-2);
    }}

    //
    //
    //
    //
    //
    // ///////////////////////////
    // // Alternate Kij calc stuff
    // ///////////////////////////
    // double LieBeta_gamma_polar[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    // double LieBeta_gamma_cart[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    // double dt_gamma_again[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    //
    // // total lie derivative
    // LieBeta_gamma_polar[0][0] = (X*betaR*dadr - 2.*a*betaR*dXdr + 2.*a*X*dbetaR_dr)/X3;
    // LieBeta_gamma_polar[1][1] = r*betaR*(r*X*dbdr + 2.*b*(X-r*dXdr))/X3;
    // LieBeta_gamma_polar[2][2] = r*sintheta2*betaR*(r*X*dbdr + 2.*b*(X-r*dXdr))/X3;
    //
    // // the beta deriv of metric (part 1) ONLY for debug
    // // LieBeta_gamma_polar[0][0] = (X*betaR*dadr - 2.*a*betaR*dXdr)/X3;
    // // LieBeta_gamma_polar[1][1] = r*betaR*(r*X*dbdr + 2.*b*(X-r*dXdr))/X3;
    // // LieBeta_gamma_polar[2][2] = r*betaR*(r*X*dbdr + 2.*b*(X-r*dXdr))/X3;
    //
    // // the deriv of beta term (part 2) ONLY for debug
    // // LieBeta_gamma_polar[0][0] = 2.*a*dbetaR_dr/X2;
    //
    // for (int i=0; i<3; i++){
    // for (int j=0; j<3; j++){
    // for (int m=0; m<3; m++){
    // for (int n=0; n<3; n++){
    //    LieBeta_gamma_cart[i][j] += dxp_dxc[n][j] * dxp_dxc[m][i] * LieBeta_gamma_polar[m][n];
    // }}}}
    //
    // for (int i=0; i<3; i++){
    // for (int j=0; j<3; j++){
    //   dt_gamma_again[i][j] = LieBeta_gamma_cart[i][j] - 2.*alpha*Kij_rest_frame[i][j];
    // }}
    //
    // ///////////////////////////////
    // // third time trying lucky ??
    // ///////////////////////////////
    //
    // double tilde_Kij[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    //
    // // boosted d_i beta^j
    // double tilde_diBi[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    //
    //
    // // d_i gamma tilde  =  LL d_i gamma
    //
    // for (int i=0; i<3; i++){
    // for (int j=0; j<3; j++){
    // for (int r=0; r<3; r++){
    // for (int s=0; s<3; s++){
    //   tilde_Kij[i][j] += Lambda[i][r] * Lambda[j][s]
    //                    * (c_ + boosted_beta_U[0] * s_)
    //                    * dt_gamma_again[r][s];
    //   tilde_Kij[i][j] += Lambda[i][r] * Lambda[j][s]
    //                    * (s_ - boosted_beta_U[0] * c_)
    //                    * d_4metric[0][r][s];
    //   tilde_Kij[i][j] += Lambda[i][r] * Lambda[j][s]
    //                    * (- boosted_beta_U[1])
    //                    * d_4metric[1][r][s];
    // }}}}
    //
    // for (int i=0; i<3; i++){
    // for (int j=0; j<3; j++){
    //   tilde_Kij[i][j] = -1./(2. * boost_lapse) * tilde_Kij[i][j];
    // }}


    // temporary storage in A_ij
    // FOR2(i,j)
    // {
    //     vars.A[i][j] += tilde_Kij[i][j];
    // }
    // //cartoon terms
    // vars.Aww += tilde_Kij[2][2];



    ///////////////////////
    // final var storage
    ///////////////////////

    current_cell.store_vars(vars);
}




#endif /* EMSBH_READ_IMPL_HPP_ */
