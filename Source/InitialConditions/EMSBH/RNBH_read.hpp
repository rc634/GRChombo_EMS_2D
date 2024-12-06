/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef RNBH_READ_HPP_
#define RNBH_READ_HPP_


#include "RNBHSolution.hpp"
#include "Cell.hpp"
#include "EMSCouplingFunction.hpp"
// #include "EinsteinMaxwellDilatonField.hpp" // not sure if needed
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "parstream.H" //gives pout
#include "simd.hpp"
#include <vector>

//! Class which solves for the initial data for a spherically symmetric emdbh

class RNBH_read
{

  public:
    //! The constructor
    RNBH_read(EMSBH_params_t a_params_EMSBH,
              CouplingFunction::params_t a_params_coupling_function,
              double a_G_Newton, double a_dx, int a_verbosity);

    //! Computes the 1d solution and stores in m_1d_sol
    void compute_1d_solution();

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;



  protected:

    double m_dx;
    double m_G_Newton;
    EMSBH_params_t m_params_EMSBH;  //!< The complex scalar field params
    CouplingFunction::params_t m_params_coupling_function; //!< The potential params
    int m_verbosity;

  public:

    RNBHSolution m_1d_sol; /*<
    The object that stores the solution found by reading data file */


};

#include "RNBH_read.impl.hpp"

#endif /* EMSBH_READ_HPP_ */
