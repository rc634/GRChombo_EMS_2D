/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(INITIALSCALARDATA_2D_HPP_)
#error "This file should only be included through InitialScalarData_2D.hpp"
#endif

#ifndef INITIALSCALARDATA_2D_IMPL_HPP_
#define INITIALSCALARDATA_2D_IMPL_HPP_

// Compute the value of the initial vars on the grid
template <class data_t>
void InitialScalarData_2D::compute(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below

    // Get coords and radius
    Coordinates<data_t> coords(current_cell, m_dx, m_init_SF_params.centerSF);
    data_t rr = coords.get_radius();

    // data_t phi = m_init_SF_params.phi0 +
    //              m_init_SF_params.amplitude *
    //                  (cos(2. * M_PI * coords.x / m_init_SF_params.L) +
    //                   cos(2. * M_PI * coords.y / m_init_SF_params.L));
    data_t width = m_init_SF_params.L;
    data_t phi = m_init_SF_params.amplitude * exp(-rr * rr / (width * width)) +
                 m_init_SF_params.amplitude * .1;
    if (!(m_init_SF_params.centerSF_2[0] < 0))
    {
        Coordinates<data_t> coords_2(current_cell, m_dx,
                                     m_init_SF_params.centerSF_2);
        data_t rr_2 = coords_2.get_radius();
        phi +=
            m_init_SF_params.amplitude * exp(-rr_2 * rr_2 / (width * width)) +
            m_init_SF_params.amplitude * .1;
    }

    data_t Pi = 0;

    vars.phi = phi;
    vars.Pi = Pi;
    vars.lapse = 1.;
    vars.chi = 1.;
    FOR(i) { vars.h[i][i] = 1.; }
    vars.hww = 1.;
    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* INITIALSCALARDATA_2D_IMPL_HPP_ */
