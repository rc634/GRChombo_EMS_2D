// Gullstrand-Painlevel coordinates
#if !defined(BLACKSTRING_HPP_)
#error "This file should only be included through BlackString.hpp"
#endif

#ifndef BLACKSTRING_IMPL_HPP_
#define BLACKSTRING_IMPL_HPP_

#include "BlackString.hpp"
#include "CCZ4CartoonVars.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

//#include <math.h>

// Computes Black String Solution
template <class data_t>
void BlackString::compute(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below
    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

    // work out where we are on the grid
    const data_t x = coords.x;
    const double cartoon_coord = coords.y;
    const double Lhor = m_params.Lhor;
    const double regFrac = m_params.regFrac;
    const double pert_amp = m_params.pert_amp;
    const double pert_freq = m_params.pert_freq;
    const int L = m_params.L; // could be a double?
    const int N1 = m_params.N1;
    const int N2 = m_params.N2;
    const double rHor = Lhor;
    // turducken
    auto cartoon_coord_regularised = cartoon_coord;
    if (cartoon_coord < regFrac * rHor)
    {
        cartoon_coord_regularised = regFrac * rHor;
    }

    //    data_t chi = 1.0 + pert_amp * cos(x*2*M_PI*N2*pert_freq/(L*N1))*
    //    exp(-pow(cartoon_coord/rHor-rHor/cartoon_coord,2)); data_t chi = 1.0 +
    //    pert_amp * exp(-pert_freq*pow(N1*L/N2/2.0-x,2)) *
    //    exp(-pow(cartoon_coord/rHor-rHor/cartoon_coord,2));
    data_t chi =
        1.0 + pert_amp * sin(x * 2 * M_PI * N2 * pert_freq / (L * N1)) *
                  exp(-pow(cartoon_coord / rHor - rHor / cartoon_coord, 2));

    vars.chi = chi;

    vars.h[0][0] = 1.0;

    vars.h[1][1] = 1.0;

    vars.hww = 1.0;

    vars.K = 3.0 * pow(rHor / cartoon_coord_regularised, 0.5) /
             (2.0 * cartoon_coord_regularised);

    vars.A[0][0] = -3.0 / 4.0 * pow(rHor / cartoon_coord_regularised, 0.5) /
                   (2.0 * cartoon_coord_regularised);

    vars.A[1][1] = -7.0 / 4.0 * pow(rHor / cartoon_coord_regularised, 0.5) /
                   (2.0 * cartoon_coord_regularised);

    vars.Aww = 5.0 / 4.0 * pow(rHor / cartoon_coord_regularised, 0.5) /
               (2.0 * cartoon_coord_regularised);

    //    vars.Gamma[1] = 0.0; //By a separate class/compute function

    // Populate the variables on the grid
    // NB We stil need to set Gamma^i which is NON ZERO
    // but we do this via a separate class/compute function
    // as we need the gradients of the metric which are not yet available

    // lapse
    switch (m_initial_lapse)
    {
    case 0:
        vars.lapse = 1.0;
        break;
    case 1:
        vars.lapse = exp(-1.0 / cartoon_coord_regularised);
        break;
    default:
        MayDay::Error("BlackString::Supplied initial lapse not supported.");
    }
    // Shift
    switch (m_initial_shift)
    {
    case 0:
        break;
    case 1:
        vars.shift[1] = 1.0 * pow(rHor / cartoon_coord_regularised, 0.5);
        break;
    default:
        MayDay::Error("BlackString::Supplied initial shift not supported.");
    }
    //    vars.B[1]=1.0/1.5/pow(y_regularised,2.);

    current_cell.store_vars(vars);
}

#endif /* BLACKSTRING_IMPL_HPP_ */
