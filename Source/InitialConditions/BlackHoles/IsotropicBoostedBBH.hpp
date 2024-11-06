/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ISOTROPICBOOSTEDBBH_HPP_
#define ISOTROPICBOOSTEDBBH_HPP_

#include "IsotropicBoostedBH.hpp"

// Superposition of 2 IsotropicBoostedBH solutions
// For simplicity, it only boosts the BH in the 'x' direction of the momentum
// (so only bh_params.momentum[0] matters)
// Note: the momentum parameter of BoostedBH::params_t in this class is used as
// velocity, not momentum, and should be between ]-1,1[
class IsotropicBoostedBBH
{
  public:
    //! The constructor
    IsotropicBoostedBBH(BoostedBH::params_t a_bh1_params,
                        BoostedBH::params_t a_bh2_params, double a_dx)
        : m_bh1(a_bh1_params, m_dx), m_bh2(a_bh2_params, m_dx), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    IsotropicBoostedBH m_bh1;
    IsotropicBoostedBH m_bh2;
    double m_dx;
};

#include "IsotropicBoostedBBH.impl.hpp"

#endif /* ISOTROPICBOOSTEDBBH_HPP_ */
