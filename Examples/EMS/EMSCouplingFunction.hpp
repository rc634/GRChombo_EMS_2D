/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COUPLINGFUNCTION_HPP_
#define COUPLINGFUNCTION_HPP_

//#include "simd.hpp"

// coupilng function of lagangean, defined by
// ~e^{-2 alpha f(phi)} F^2
// with f = f0 + f1 * phi + f2 * phi * phi


class CouplingFunction
{
  public:
    struct params_t
    {
        double alpha;
        double f0;
        double f1;
        double f2;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    CouplingFunction(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the phi field here
    template <class data_t, template <typename> class vars_t>
    void compute_coupling(data_t &f_of_phi,
                           data_t &f_prime_of_phi,
                           data_t &coupling_of_phi,
                           const vars_t<data_t> &vars) const
    {
        // calculate f(phi) and then teh exponential coupling
        f_of_phi = m_params.f0 + m_params.f1 * vars.phi
                               + m_params.f2 * vars.phi * vars.phi;
        f_prime_of_phi = m_params.f1 + 2. * m_params.f2 * vars.phi;

        f_of_phi *= m_params.alpha;
        f_prime_of_phi *= m_params.alpha;

        coupling_of_phi = exp(-2.*f_of_phi);

    }


};



// standalone function
template <class data_t>
void compute_coupling_of_phi(
                       data_t a_alpha,
                       data_t a_f0,
                       data_t a_f1,
                       data_t a_f2,
                       data_t &f_of_phi,
                       data_t &f_prime_of_phi,
                       data_t &coupling_of_phi,
                       data_t &a_phi)
{
    // calculate f(phi) and then teh exponential coupling
    f_of_phi = a_f0 + a_f1 * a_phi
                           + a_f2 * a_phi * a_phi;
    f_prime_of_phi = a_f1 + 2. * a_f2 * a_phi;

    f_of_phi *= a_alpha;
    f_prime_of_phi *= a_alpha;

    coupling_of_phi = exp(-2.*f_of_phi);

}

#endif /* COUPLINGFUNCTION_HPP_ */
