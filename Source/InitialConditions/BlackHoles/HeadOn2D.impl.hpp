// Gullstrand-Painlevel coordinates
#if !defined(BLACKSTRING_HPP_)
#error "This file should only be included through HeadOn2D.hpp"
#endif

#ifndef BLACKSTRING_IMPL_HPP_
#define BLACKSTRING_IMPL_HPP_

#include "CCZ4CartoonVars.hpp"
#include "HeadOn2D.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

//#include <math.h>

// Computes Black String Solution
template <class data_t> void HeadOn2D::compute(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below
    Coordinates<data_t> coords1(current_cell, m_dx, m_bh1_params.center);
    Coordinates<data_t> coords2(current_cell, m_dx, m_bh2_params.center);

    // work out where we are on the grid
    const double m1 = m_bh1_params.mass;
    const data_t r1 = coords1.get_radius();
    const double m2 = m_bh2_params.mass;
    const data_t r2 = coords2.get_radius();

    //    data_t chi = 1.0 + pert_amp * cos(x*2*M_PI*N2*pert_freq/(L*N1))*
    //    exp(-pow(z/rHor-rHor/z,2)); data_t chi = 1.0 + pert_amp *
    //    exp(-pert_freq*pow(N1*L/N2/2.0-x,2)) * exp(-pow(z/rHor-rHor/z,2));
    data_t chi = 1.0 + m1 / (2 * r1) + m2 / (2 * r2);

    vars.chi = pow(chi, -4);

    vars.h[0][0] = 1.0;

    vars.h[1][1] = 1.0;

    vars.hww = 1.0;

    vars.K = 0;

    vars.A[0][0] = 0;

    vars.A[1][1] = 0;

    vars.Aww = 0;

    //    vars.Gamma[1] = 0.0; //By a separate class/compute function

    // Populate the variables on the grid
    // NB We stil need to set Gamma^i which is NON ZERO
    // but we do this via a separate class/compute function
    // as we need the gradients of the metric which are not yet available

    // lapse
    vars.lapse = 1.0;
    current_cell.store_vars(vars);
}

#endif /* BLACKSTRING_IMPL_HPP_ */
