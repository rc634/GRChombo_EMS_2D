/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_2D_HPP_
#define INITIALSCALARDATA_2D_HPP_

#include "CCZ4CartoonVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "parstream.H" //gives pout
#include "Potential.hpp"
#include "simd.hpp"
#include <array>
#include <vector>
#include "SimulationParameters.hpp"

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class InitialScalarData_2D
{

  public:
    struct params_t
    {
        double phi0;
        double amplitude;
        double k;
        double L;
        std::array<double, CH_SPACEDIM>
            centerSF;
        std::array<double, CH_SPACEDIM>
            centerSF_2;
    };
    //! The constructor
    InitialScalarData_2D(params_t a_init_SF_params, Potential::params_t a_potential_params, double a_dx)
        : m_init_SF_params(a_init_SF_params), m_potential_params(a_potential_params), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const params_t m_init_SF_params;
    const Potential::params_t m_potential_params;
    double m_dx;
};

#include "InitialScalarData_2D.impl.hpp"

#endif /* INITIALSCALARDATA_2D_HPP_ */
