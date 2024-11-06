/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ISOTROPICBOOSTEDBH_HPP_
#define ISOTROPICBOOSTEDBH_HPP_

#include "Cell.hpp"

// This is a simple class from the 2D code that was useful to test the AHFinder
// for high boosts, as it produces highly deformed (ellipsoidal) AHs
// For simplicity, it only boosts the BH in the 'x' direction of the momentum
// (so only bh_params.momentum[0] matters)
// Note: the momentum parameter of BoostedBH::params_t in this class is used as
// velocity, not momentum, and should be between ]-1,1[
class IsotropicBoostedBH
{
  public:
#if CH_SPACEDIM == 3
    template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;
#elif CH_SPACEDIM == 2
    template <class data_t> using Vars = CCZ4CartoonVars::VarsWithGauge<data_t>;
#endif

    //! The constructor
    IsotropicBoostedBH(BoostedBH::params_t a_bh_params, double a_dx)
        : m_bh_params(a_bh_params), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    template <class data_t>
    void computeMetric(Tensor<2, data_t, 3> &gammaLL,
                       Tensor<2, data_t, 3> &gammaUU, data_t &shift_x,
                       Tensor<2, data_t, 3> &KLL,
                       const Coordinates<data_t> &coords) const;

    template <class data_t>
    static void computeVars(Vars<data_t> &vars,
                            const Tensor<2, data_t, 3> &gammaLL,
                            const Tensor<2, data_t, 3> &gammaUU,
                            data_t &shift_x, const Tensor<2, data_t, 3> &KLL);

    const BoostedBH::params_t m_bh_params;

  protected:
    double m_dx;
};

#include "IsotropicBoostedBH.impl.hpp"

#endif /* ISOTROPICBOOSTEDBH_HPP_ */
