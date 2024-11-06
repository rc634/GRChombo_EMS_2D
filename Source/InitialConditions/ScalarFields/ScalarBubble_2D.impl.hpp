/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCALARBUBBLE_2D_HPP_)
#error "This file should only be included through ScalarBubble_2D.hpp"
#endif

#ifndef SCALARBUBBLE_2D_IMPL_HPP_
#define SCALARBUBBLE_2D_IMPL_HPP_

// Compute the value of the initial vars on the grid
template <class data_t>
void ScalarBubble_2D::compute(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below

    // Get coords and radius
    Coordinates<data_t> coords(current_cell, m_dx, m_bubble_params.centerSF);
    data_t rr_1 = coords.get_radius(); // distance to center of bubble 1

    // Find parameters dependent on potential coefficients
    data_t m = m_potential_params.phi0 * sqrt(m_potential_params.lambda);
    data_t radius_critical =
        2 * m_potential_params.lambda / (3 * m_potential_params.epsilon * m);
    data_t R0 = radius_critical * m_bubble_params.radius_factor;
    data_t gamma = R0 * m_potential_params.epsilon * m / m_potential_params.lambda +
                   4 * m_potential_params.lambda * m_potential_params.lambda /
                       (27. * R0 * R0 * m_potential_params.epsilon *
                        m_potential_params.epsilon * m * m);
    data_t v = sqrt((gamma * gamma - 1) / (gamma * gamma));

    data_t phi, Pi;
    if (m_bubble_params.centerSF_2[0] < 0.)
    {
        phi = compute_phi(gamma, m, rr_1, R0);
        Pi = compute_Pi(gamma, m, v, rr_1, R0);
    }
    else
    {
        // This doesn't currently work
        Coordinates<data_t> coords_2(current_cell, m_dx,
                                     m_bubble_params.centerSF_2);
        data_t rr_2 = coords_2.get_radius(); // distance to center of bubble 2
        data_t dist = simd_min(rr_1, rr_2);
        phi = compute_phi(gamma, m, dist, R0);
        Pi = compute_Pi(gamma, m, v, dist, R0);
    }

    vars.phi = phi;
    vars.Pi = Pi;
    vars.lapse = 1.;
    vars.chi = 1.;
    FOR(i) { vars.h[i][i] = 1.; }
    vars.hww = 1.;
    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

template <class data_t>
data_t ScalarBubble_2D::compute_phi(data_t gamma, data_t m, data_t rr,
                                    data_t R0) const
{
    return (m_potential_params.false_vacuum - m_potential_params.true_vacuum) *
               tanh(gamma * m * (rr - R0) / 2.) / 2. +
           (m_potential_params.false_vacuum + m_potential_params.true_vacuum) / 2.;
}

template <class data_t>
data_t ScalarBubble_2D::compute_Pi(data_t gamma, data_t m, data_t v, data_t rr,
                                   data_t R0) const
{
    return -m_potential_params.phi0 * gamma * m * v *
           (1 - tanh(gamma * m * (rr - R0) / 2.) *
                    tanh(gamma * m * (rr - R0) / 2.)) /
           2.;
}

#endif /* SCALARBUBBLE_2D_IMPL_HPP_ */
