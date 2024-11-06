/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMDEXTRACTIONTAGGINGCRITERION_HPP_
#define EMDEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "EinsteinMaxwellScalarField.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SimulationParametersBase.hpp"
#include "Tensor.hpp"

class EMDExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const extraction_params_t m_params;
    const int m_level;
    const double m_threshold_A;
    const double m_threshold_phi;
    const double m_threshold_chi;

    template <class data_t>
    using MatterVars =
                  typename EinsteinMaxwellScalarField<>::template Vars<data_t>;

    // Vars object for chi
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

  public:
    EMDExtractionTaggingCriterion(
        const double a_dx, const int a_level,
        const extraction_params_t a_params, const double a_threshold_A,
        const double a_threshold_phi, const double a_threshold_chi)
        : m_dx(a_dx), m_deriv(a_dx), m_params(a_params), m_level(a_level),
                                         m_threshold_A(a_threshold_A),
                                         m_threshold_phi(a_threshold_phi),
                                         m_threshold_chi(a_threshold_chi){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);
        const auto d1chi = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<MatterVars>(current_cell);
        const auto d2chi = m_deriv.template diff2<Vars>(current_cell);

        data_t d2_phi_ratio = 0.;
        data_t d2_chi_ratio = 0.;
        data_t d2_A_ratio = 0.;

        FOR2(i, j)
        {
            //data_t mod_d2_A_ij = d2.Xi[i][j] * d2.Xi[i][j];
            data_t mod_d2_Pi_ij = d2.Pi[i][j] * d2.Pi[i][j];
            data_t mod_d2_phi_ij = d2.phi[i][j] * d2.phi[i][j];

            //data_t abs_d1d1_A_ij = abs(d1.Xi[i] * d1.Xi[j]);
            data_t abs_d1d1_Pi_ij = abs(d1.Pi[i] * d1.Pi[j]);
            data_t abs_d1d1_phi_ij = abs(d1.phi[i] * d1.phi[j]);

            //d2_A_ratio += mod_d2_A_ij / (abs_d1d1_A_ij + 1e-5);
            d2_phi_ratio += mod_d2_Pi_ij / (abs_d1d1_Pi_ij + 1e-5) +
                            mod_d2_phi_ij / (abs_d1d1_phi_ij + 1e-5);
            d2_chi_ratio += d2chi.chi[i][j] * d2chi.chi[i][j] /
                            (1e-2 + abs(d1chi.chi[i] * d1chi.chi[j]));
        }
        data_t hack_coeff = 1.; // mega hacked
        data_t numerator = 0.3*pow(m_level,2) + 0.1*pow(m_level,3);
        data_t denomenator = 1. + 0.000006*pow(m_level,7);
        data_t hack_dx = m_dx * pow(2,m_level)/(pow(2,m_level)-numerator/denomenator);

        data_t criterion_phi = hack_dx * sqrt(d2_phi_ratio) / m_threshold_phi;
        data_t criterion_chi = hack_dx * sqrt(d2_chi_ratio) / m_threshold_chi;
        //data_t criterion_A = hack_dx * sqrt(d2_A_ratio) / m_threshold_A;

        data_t criterion = simd_max(criterion_chi, criterion_phi);
        //criterion = simd_max(criterion, criterion_A);

        for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
        {
            // regrid if within extraction level and not at required refinement
            if (m_level < m_params.extraction_levels[iradius])
            {
                const Coordinates<data_t> coords(current_cell, m_dx,
                                                 m_params.extraction_center);
                const data_t r = coords.get_radius();
                // add a 20% buffer to extraction zone so not too near to
                // boundary
                auto regrid = simd_compare_lt(
                    r, 1.2 * m_params.extraction_radii[iradius]);
                criterion = simd_conditional(regrid, 100.0, criterion);
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* EMDEXTRACTIONTAGGINGCRITERION_HPP_ */
