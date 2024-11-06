/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef LAYERTAGGINGCRITERION_HPP_
#define LAYERTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"

class LayerTaggingCriterion
{
  protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_y_regrid_lim;

    // const FourthOrderDerivatives m_deriv;
    // const double m_threshold_phi;
    // const double m_threshold_K;
    // const double m_threshold_rho;
    // const double m_L;
    // const int m_level;
    // const double m_max_rho_radius;
    // const double m_regrid_wall_width;
    // const int m_wall_min_regrid_level;
    // const int m_away_max_regrid_level;
    // const double m_time;

  public:
    LayerTaggingCriterion(double dx, std::array<double, CH_SPACEDIM> a_center,
                           double a_y_regrid_lim
                           // , double threshold_phi, double threshold_K,
                           // double threshold_rho,
                           // double a_L, int a_level, double a_max_rho_radius,
                           // double a_regrid_wall_width,
                           // double a_wall_min_regrid_level,
                           // double a_away_max_regrid_level,
                           // double time
                           )
        : m_dx(dx), m_center(a_center), m_y_regrid_lim(a_y_regrid_lim)
        //   ,
        //   m_threshold_K(threshold_K), m_threshold_rho(threshold_rho), m_L(a_L),
        //   m_level(a_level), m_max_rho_radius(a_max_rho_radius),
        //   m_regrid_wall_width(a_regrid_wall_width),
        //   m_wall_min_regrid_level(a_wall_min_regrid_level),
        //   m_away_max_regrid_level(a_away_max_regrid_level), m_time(time)
          {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        data_t criterion = 0.0;

        if (coords.y < m_y_regrid_lim)
        {
            criterion = 10000.0;
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* LAYERTAGGINGCRITERION_HPP_ */
