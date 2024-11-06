/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEFAULTMEDCOUPLINGFUNCTION_HPP_
#define DEFAULTMEDCOUPLINGFUNCTION_HPP_

#include "Tensor.hpp"
#include "simd.hpp"

class DefaultEMDCouplingFunction
{
  public:
    //! The constructor
    DefaultEMDCouplingFunction() {}

    //! Set the potential function for the emd field here to unity
    template <class data_t, template <typename> class vars_t>
    void compute_coupling(data_t &f_of_phi,
                           data_t &f_prime_of_phi,
                           data_t &coupling_of_phi,
                           const vars_t<data_t> &vars) const
    {
        f_of_phi=0.0;
        f_prime_of_phi=0.0;
        coupling_of_phi=1.0;
    }
};

#endif /* DEFAULTMEDCOUPLINGFUNCTION_HPP_ */
