/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXPOTENTIAL_HPP_
#define COMPLEXPOTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
public:
    struct params_t
    {
        // double scalar_mass;
        double lambda;
        double epsilon;
        double phi0;
        double Vbasetrue;
        double true_vacuum;
        double false_vacuum;

        // Lewicki params
        bool use_Lewicki_init_data;
        double L_pot_a;
        double L_pot_v;

        //From Robin's code
        double scalar_mass;
        double phi4_coeff;
        bool solitonic;
        double sigma_soliton;

    };

private:
    params_t m_params;

public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}



    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_modulus_phi_squared,
                           data_t &dVdmodulus_phi_squared,
                           const vars_t<data_t> &vars) const
    {
        // First calculate |phi|^2
        data_t modulus_phi_squared =
                vars.phi * vars.phi + vars.phi_Im * vars.phi_Im;

        if (!m_params
                .solitonic) // if the star is not solitonic, make lambda star
        {
            // The potential value at phi (note the convention with factors of
            // 1/2) m^2 |phi|^2 + lambda/2 |phi|^4
            V_of_modulus_phi_squared =
                    m_params.scalar_mass * m_params.scalar_mass *
                    modulus_phi_squared +
                    0.5 * m_params.phi4_coeff * modulus_phi_squared *
                    modulus_phi_squared;

            //   V_of_phi = V_of_modulus_phi_squared; maybe we dont need this

            // The potential gradient at phi
            // m^2 + lambda |phi|^2
            dVdmodulus_phi_squared =
                    m_params.scalar_mass * m_params.scalar_mass +
                    m_params.phi4_coeff * modulus_phi_squared;

            //   dVdphi = dVdmodulus_phi_squared; maybe we dont need this
        }
        else // else star is solitonic
        {
            // The potential value at phi (note the convention with factors of
            // 1/2) m^2 |phi|^2 * (1 - \frac{2|phi|^2}{sigma^2})^2
            V_of_modulus_phi_squared =
                    m_params.scalar_mass * m_params.scalar_mass *
                    modulus_phi_squared *
                    pow(1. - 2. * modulus_phi_squared /
                             (m_params.sigma_soliton * m_params.sigma_soliton),
                        2);

            //   V_of_phi = V_of_modulus_phi_squared; maybe we dont need this

            dVdmodulus_phi_squared =
                    m_params.scalar_mass * m_params.scalar_mass *
                    (1. - 6. * modulus_phi_squared /
                          (m_params.sigma_soliton * m_params.sigma_soliton)) *
                    (1. - 2. * modulus_phi_squared /
                          (m_params.sigma_soliton * m_params.sigma_soliton));

            //   dVdphi = dVdmodulus_phi_squared; maybe we dont need this
        }
    }
};

#endif /* COMPLEXPOTENTIAL_HPP_ */

