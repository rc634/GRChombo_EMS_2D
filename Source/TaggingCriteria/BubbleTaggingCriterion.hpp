/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BUBBLETAGGINGCRITERION_HPP_
#define BUBBLETAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "ComplexScalarField.hpp"
#include "Coordinates.hpp"
#include "SimulationParametersBase.hpp"

class BubbleTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_phi;
    const double m_threshold_chi;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_L;
    const int m_level;
    const double m_max_rho_radius;
    const double m_regrid_wall_width;
    const int m_wall_min_regrid_level;
    const int m_away_max_regrid_level;
    const double m_time;
    const extraction_params_t m_params;


    template <class data_t>
    using MatterVars = typename ComplexScalarField<>::template Vars<data_t>;

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
    BubbleTaggingCriterion(double dx, double a_threshold_phi, double a_threshold_chi,
                            std::array<double, CH_SPACEDIM> a_center,
                            double a_L, int a_level, double a_max_rho_radius,
                            double a_regrid_wall_width,
                            double a_wall_min_regrid_level,
                            double a_away_max_regrid_level,
                            double time,
                            const extraction_params_t a_params)
        : m_dx(dx), m_deriv(dx), m_threshold_phi(a_threshold_phi),
          m_threshold_chi(a_threshold_chi), m_center(a_center), m_L(a_L),
          m_level(a_level), m_max_rho_radius(a_max_rho_radius),
          m_regrid_wall_width(a_regrid_wall_width),
          m_wall_min_regrid_level(a_wall_min_regrid_level),
          m_away_max_regrid_level(a_away_max_regrid_level), m_time(time), m_params(a_params){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);
        const auto d1chi = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<MatterVars>(current_cell);
        const auto d2chi = m_deriv.template diff2<Vars>(current_cell);

        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        
        // Standard criterion definitions
 //       data_t mod_d1_rho = 0.0; 
        data_t d2_phi_ratio = 0.;
        data_t d2_chi_ratio = 0.;

        FOR2(idir, jdir)
        {
            data_t mod_d2_Pi_ij = d2.Pi[idir][jdir] * d2.Pi[idir][jdir] +
                                  d2.Pi_Im[idir][jdir] * d2.Pi_Im[idir][jdir];
            data_t mod_d2_phi_ij =
                    d2.phi[idir][jdir] * d2.phi[idir][jdir] +
                    d2.phi_Im[idir][jdir] * d2.phi_Im[idir][jdir];
            data_t abs_d1d1_Pi_ij = abs(d1.Pi[idir] * d1.Pi[jdir] +
                                        d1.Pi_Im[idir] * d1.Pi_Im[jdir]);
            data_t abs_d1d1_phi_ij = abs(d1.phi[idir] * d1.phi[jdir] +
                                         d1.phi_Im[idir] * d1.phi_Im[jdir]);
            d2_phi_ratio += mod_d2_Pi_ij / (abs_d1d1_Pi_ij + 1e-5) +
                            mod_d2_phi_ij / (abs_d1d1_phi_ij + 1e-5);
            d2_chi_ratio += d2chi.chi[idir][jdir] * d2chi.chi[idir][jdir] /
                            (1e-2 + abs(d1chi.chi[idir] * d1chi.chi[jdir]));
        }

        data_t criterion_phi = m_dx * sqrt(d2_phi_ratio) / m_threshold_phi;
        data_t criterion_chi = m_dx * sqrt(d2_chi_ratio) / m_threshold_chi;

        data_t criterion = simd_max(criterion_chi, criterion_phi);

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

 /*
        // Hacking the shit out of this: can't regrid wrt rho at t=0, so using
        // phi instead
        data_t mod_d1_phi = 0.0;
        // if (m_time == 0.0)
        // {
            Tensor<1, data_t> d1_phi;
            FOR(idir)
            {
                m_deriv.diff1(d1_phi, current_cell, idir, c_phi);
                mod_d1_phi += d1_phi[idir] * d1_phi[idir];
            }
        // }


        data_t criterion = m_dx * (sqrt(mod_d1_rho) / m_threshold_rho +
                                   sqrt(mod_d1_K) / m_threshold_K +
                                   sqrt(mod_d1_phi) / m_threshold_phi);

*/

        // if (rxy > min(m_L / 2.0, 1.2 * m_max_rho_radius) && m_level > m_away_max_regrid_level - 1)
        // {
        //     criterion = 0.0;
        // }
        
        // Making sure we regrid everywhere within a distance
        // m_regrid_wall_width from the wall. The code thinks the wall is at
        // r_max, where rho is maximum. Therefore, we only want to do this when
        // r_max < L.
        // if (simd_compare_lt(m_max_rho_radius - m_regrid_wall_width, rxy) &&
        //     simd_compare_lt(rxy, m_max_rho_radius + m_regrid_wall_width) &&
        //     m_level < m_wall_min_regrid_level - 0 && m_max_rho_radius < m_L)
        // {
        //     // pout() << "setting criterion to 1 at radius " << rxy  << " at
        //     // level " << m_level << endl;
        //     criterion = 10000.0;
        // }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* BUBBLETAGGINGCRITERION_HPP_ */
