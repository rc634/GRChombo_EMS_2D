/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARBUBBLE_2D_HPP_
#define SCALARBUBBLE_2D_HPP_

#include "CCZ4CartoonVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "SimulationParameters.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "parstream.H" //gives pout
#include "Potential.hpp"
#include "simd.hpp"
#include <array>
#include <vector>

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class ScalarBubble_2D
{

  public:
    struct params_t
    {
        double radius_factor;
        std::array<double, CH_SPACEDIM> centerSF;
        std::array<double, CH_SPACEDIM> centerSF_2;
    };
    //! The constructor
    ScalarBubble_2D(params_t a_bubble_params, Potential::params_t a_potential_params, double a_dx)
        : m_bubble_params(a_bubble_params), m_potential_params(a_potential_params), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    template <class data_t>
    data_t compute_phi(data_t gamma, data_t m, data_t rr, data_t R0) const;
    template <class data_t>
    data_t compute_Pi(data_t gamma, data_t m, data_t v, data_t rr, data_t R0) const;

  protected:
    const params_t m_bubble_params;
    const Potential::params_t m_potential_params;
    double m_dx;
};

#include "ScalarBubble_2D.impl.hpp"

#endif /* SCALARBUBBLE_2D_HPP_ */
