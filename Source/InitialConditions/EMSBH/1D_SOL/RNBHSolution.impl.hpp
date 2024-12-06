/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(RNBHSOLUTION_HPP_)
#error "This file should only be included through RNBHSolution.hpp"
#endif

#ifndef RNBHSOLUTION_IMPL_HPP_
#define RNBHSOLUTION_IMPL_HPP_

RNBHSolution::RNBHSolution() {}


double RNBHSolution::calc_Er(const double r) const
{
    // create the Reissner Nordstrom solution
    double out=0.;
    double rep = sqrt(8.*M_PI); // REP = Root Eight Pi

    out = - Q / (rep * calc_psi(r) * r * r);

    return out;
}

double RNBHSolution::calc_psi(const double r) const
{
    // create the Reissner Nordstrom solution
    double out=0.;

    out = 1. + M / r + (M * M - Q * Q)/(4. * r * r);

    return out;
}





void RNBHSolution::set_initialcondition_params(
    EMSBH_params_t m_params_EMSBH,
    CouplingFunction::params_t m_params_coupling_function)
{     // conformal factor
    G = m_params_EMSBH.Newtons_constant;
    M = m_params_EMSBH.bh_mass;
    Q = m_params_EMSBH.bh_charge;
}

#endif /* RNBHSOLUTION_IMPL_HPP_ */
