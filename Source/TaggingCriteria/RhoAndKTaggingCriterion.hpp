/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef RHOANDKTAGGINGCRITERION_HPP_
#define RHOANDKTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class RhoAndKTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_phi;
    const double m_threshold_K;
    const double m_threshold_rho;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_L;
    const int m_level;
    const double m_max_rho_radius;
    const double m_regrid_wall_width;
    const int m_wall_min_regrid_level;
    const int m_away_max_regrid_level;
    const double m_time;

  public:
    RhoAndKTaggingCriterion(double dx, double threshold_phi, double threshold_K,
                            double threshold_rho,
                            std::array<double, CH_SPACEDIM> a_center,
                            double a_L, int a_level, double a_max_rho_radius,
                            double a_regrid_wall_width,
                            double a_wall_min_regrid_level,
                            double a_away_max_regrid_level,
                            double time)

        : m_dx(dx), m_deriv(dx), m_threshold_phi(threshold_phi),
          m_threshold_K(threshold_K), m_threshold_rho(threshold_rho),
          m_center(a_center), m_L(a_L), m_level(a_level),
          m_max_rho_radius(a_max_rho_radius),
          m_regrid_wall_width(a_regrid_wall_width),
          m_wall_min_regrid_level(a_wall_min_regrid_level),
          m_away_max_regrid_level(a_away_max_regrid_level), m_time(time){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Tensor<1, data_t> d1_rho;
        FOR(idir) m_deriv.diff1(d1_rho, current_cell, idir, c_rho);

        Tensor<1, data_t> d1_K;
        FOR(idir) m_deriv.diff1(d1_K, current_cell, idir, c_K);

        data_t mod_d1_rho = 0;
        data_t mod_d1_K = 0;
        FOR(idir)
        {

            mod_d1_rho += d1_rho[idir] * d1_rho[idir];
            mod_d1_K += d1_K[idir] * d1_K[idir];
        }

        // Hacking the shit out of this: can't regrid wrt rho at t=0, so using
        // phi instead
        data_t mod_d1_phi = 0;
        if (m_time == 0.)
        {
            Tensor<1, data_t> d1_phi;
            FOR(idir)
            {
                m_deriv.diff1(d1_phi, current_cell, idir, c_phi);
                mod_d1_phi += d1_phi[idir] * d1_phi[idir];
            }
        }

        data_t criterion = m_dx * (sqrt(mod_d1_rho) / m_threshold_rho +
                                   sqrt(mod_d1_K) / m_threshold_K +
                                   sqrt(mod_d1_phi) / m_threshold_phi);

        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        data_t rxy =
            sqrt(coords.x * coords.x + coords.y * coords.y + .00000001);
        // For boson stars, don't regrid away from the star
        if (simd_compare_gt(rxy, 15.) && m_level > m_away_max_regrid_level - 1)
        {
            // pout() << "setting criterion to 1 at radius " << rxy  << " at
            // level " << m_level << endl;
            criterion = 0.0;
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* RHOANDKTAGGINGCRITERION_HPP_ */
