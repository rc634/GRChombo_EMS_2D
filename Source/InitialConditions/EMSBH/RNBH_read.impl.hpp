/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(RNBH_READ_HPP_)
#error "This file should only be included through RNBH_read.hpp"
#endif

#ifndef RNBH_READ_IMPL_HPP_
#define RNBH_READ_IMPL_HPP_

//#include "RNBHSolution_read.hpp" //for RNBHSolution class

inline RNBH_read::RNBH_read(EMSBH_params_t a_params_EMSBH,
                    CouplingFunction::params_t a_params_coupling_function,
                            double a_G_Newton, double a_dx, int a_verbosity)
    : m_dx(a_dx), m_G_Newton(a_G_Newton),
      m_params_EMSBH(a_params_EMSBH),
      m_params_coupling_function(a_params_coupling_function), m_verbosity(a_verbosity)
{
}

void RNBH_read::compute_1d_solution()
{
    try
    {
        // set initial parameters and then run the solver
        pout() << "run m_1d_sol.main()" << endl;
        m_1d_sol.set_initialcondition_params(m_params_EMSBH, m_params_coupling_function);
        pout() << "completed m_1d_sol.main()" << endl;
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t> void RNBH_read::compute(Cell<data_t> current_cell) const
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
    double psi = m_1d_sol.calc_psi(r);
    double Er = m_1d_sol.calc_Er(r);


    // spatial metrics
    double gamma_polar[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    gamma_polar[0][0] = psi*psi;
    gamma_polar[1][1] = r*r*psi*psi;;
    gamma_polar[2][2] = gamma_polar[1][1]*pow(sintheta,2);
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
        //gamma_inv[i][j] += gamma_polar_inv[m][n]*dxc_dxp[m][i]*dxc_dxp[n][j];
        gamma_inv[i][j] += gamma_polar_inv[m][n]*dxc_dxp[i][m]*dxc_dxp[j][n];
    }}}}

    // double bh_chi = pow(1. + 0.5/r,-4);
    // set gauge vars
    vars.lapse = 1./m_1d_sol.calc_psi(r);
    // vars.lapse = sqrt(bh_chi);
    vars.shift[0] = 0.;
    vars.shift[1] = 0.;

    // set geometry
    // metric det (det_gamma)
    double det_gamma = gamma[0][0]*gamma[1][1]*gamma[2][2];
    det_gamma += 2.*gamma[0][2]*gamma[0][1]*gamma[1][2];
    det_gamma -= gamma[0][0]*pow(gamma[1][2],2);
    det_gamma -= gamma[1][1]*pow(gamma[2][0],2);
    det_gamma -= gamma[2][2]*pow(gamma[0][1],2);
    // assign real chi from power 1/3 of metric det
    vars.chi = pow(det_gamma,-1./3.);
    // vars.chi = bh_chi;

    // don't need to loop z direction
    vars.K = 0.;
    FOR2(i,j)
    {
        vars.h[i][j] = vars.chi*gamma[i][j];
        vars.A[i][j] = 0.;
    }

    //cartoon terms
    vars.hww = vars.chi*gamma[2][2];
    vars.Aww = 0.;


    // electric field
    vars.Ex = dr_dx * Er;
    vars.Ey = dr_dy * Er;
    vars.Ez = dr_dz * Er;

    if (binary)
    {
      ////////////////////////////
      // read 2nd black hole
      ////////////////////////////

      // coord objects
      x = coords.x + 0.5 * separation;
      z = 0.; // coords.z;
      y = coords.y;
      cart_coords[0] = x; // only need to reset x here

      // radii and safe (divisible) radii
      r = sqrt(x * x + y * y + z * z);
      safe_r = sqrt(x * x + y * y + z * z + 10e-20);
      rho = sqrt(x * x + y * y);
      safe_rho = sqrt(x * x + y * y + 10e-20);

      // trig functions
      sintheta = rho/safe_r;
      costheta = z/safe_r;
      sinphi = y/safe_rho;
      cosphi = x/safe_rho;

      // jacobeans
      dx_dr = cosphi*sintheta;
      dy_dr = sinphi*sintheta;
      dz_dr = costheta;
      dr_dx = (x / safe_r);
      dr_dy = (y / safe_r);
      dr_dz = (z / safe_r);

      // partial cartesian coords (i) by partial polar coords (j)
      // dxc_dxp[i][j]
      // i = {x,y,z}
      // j = {r,th,ph}
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
      psi = m_1d_sol.calc_psi(r);
      Er = m_1d_sol.calc_Er(r);


      // spatial metrics
      //gamma_polar[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      gamma_polar[0][0] = psi*psi;
      gamma_polar[1][1] = r*r*psi*psi;;
      gamma_polar[2][2] = gamma_polar[1][1]*pow(sintheta,2);
      // gamma_polar_inv[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      gamma_polar_inv[0][0] = 1./gamma_polar[0][0];
      gamma_polar_inv[1][1] = 1./gamma_polar[1][1];
      gamma_polar_inv[2][2] = 1./gamma_polar[2][2];


      // create cartesian metric from transformation of polar version
      // explicit loops here include z direction
      for (int i=0; i<3; i++){
      for (int j=0; j<3; j++){
          gamma[i][j] = 0.;
          gamma_inv[i][j] = 0.;
      for (int m=0; m<3; m++){
      for (int n=0; n<3; n++){
          gamma[i][j] += gamma_polar[m][n]*dxp_dxc[m][i]*dxp_dxc[n][j];
          //gamma_inv[i][j] += gamma_polar_inv[m][n]*dxc_dxp[m][i]*dxc_dxp[n][j];
          gamma_inv[i][j] += gamma_polar_inv[m][n]*dxc_dxp[i][m]*dxc_dxp[j][n];
      }}}}

      // double bh_chi = pow(1. + 0.5/r,-4);
      // set shift
      vars.shift[0] = 0.;
      vars.shift[1] = 0.;

      // set geometry
      // metric det (det_gamma)
      det_gamma = gamma[0][0]*gamma[1][1]*gamma[2][2];
      det_gamma += 2.*gamma[0][2]*gamma[0][1]*gamma[1][2];
      det_gamma -= gamma[0][0]*pow(gamma[1][2],2);
      det_gamma -= gamma[1][1]*pow(gamma[2][0],2);
      det_gamma -= gamma[2][2]*pow(gamma[0][1],2);
      // assign real chi from power 1/3 of metric det
      data_t chi1 = vars.chi;
      data_t chi2 = pow(det_gamma,-1./3.);
      vars.chi = (chi1*chi2)/(chi2 + chi1 - chi1*chi2);
      vars.lapse = sqrt(vars.chi);

      // don't need to loop z direction
      vars.K += 0.;
      FOR2(i,j)
      {
          vars.h[i][j] += chi2*gamma[i][j]-kroneka_delta[i][j];
          vars.A[i][j] += 0.;
      }

      //cartoon terms
      vars.hww += chi2*gamma[2][2]-1.;
      vars.Aww += 0.;

      // fix the fact that superposing probbaly broke det(h) = 1
      data_t det_h = vars.h[0][0]*vars.h[1][1]*vars.hww
                   - vars.h[0][1]*vars.h[1][0]*vars.hww;
      data_t chi_h = pow(det_h,-1./3.);

      vars.chi *= chi_h;
      vars.hww /= chi_h;
      FOR2(i,j) vars.h[i][j] /= chi_h;


      // electric field
      vars.Ex += dr_dx * Er;
      vars.Ey += dr_dy * Er;
      vars.Ez += dr_dz * Er;
    }

    // scalar field
    vars.phi = 0.;
    vars.Pi = 0.;



    ////////////////////////////
    // zeros (magnetics field and electromag constrinats)
    ////////////////////////////

    vars.Xi = 0.; // magnetic constraint var
    vars.Lambda = 0.; // leccy constraint var

    // magnetic fields
    vars.Bx = 0.;
    vars.By = 0.;
    vars.Bz = 0.;


    // Store the initial values of the variables
    current_cell.store_vars(vars);

}

#endif /* RNBH_READ_IMPL_HPP_ */
