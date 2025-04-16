/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
// Problem specific includes (robin):
#include "EMSBHParams.hpp"
#include "EMSCouplingFunction.hpp"

#ifdef USE_AHFINDER
#include "AHInitialGuess.hpp"
#endif

class SimulationParameters : public SimulationParametersBase
{
public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
        check_params();
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_threshold_A", regrid_threshold_A);
        pp.load("regrid_threshold_phi", regrid_threshold_phi);
        pp.load("regrid_threshold_chi", regrid_threshold_chi);

        pp.load("G_Newton", m_G_Newton);

        // EMSBH initial data params
        pp.load("ems_data_path", emsbh_params.data_path); // the default should fail
        pp.load("gridpoints", emsbh_params.gridpoints, 40000);
        pp.load("star_centre", emsbh_params.star_centre,
                {0.5 * L, 0.5 * L});
        pp.load("binary", emsbh_params.binary, false);
        pp.load("separation", emsbh_params.separation, 0.0);
        pp.load("bh_charge", emsbh_params.bh_charge, 0.0);
        pp.load("bh_mass", emsbh_params.bh_mass, 1.0);
        pp.load("G_Newton", emsbh_params.Newtons_constant, 1.0);

        // Coupling params
        pp.load("ems_alpha", coupling_function_params.alpha, 0.0);
        pp.load("ems_f0", coupling_function_params.f0, 0.0);
        pp.load("ems_f1", coupling_function_params.f1, 0.0);
        pp.load("ems_f2", coupling_function_params.f2, 0.0);

        // Perturbation shell params
        pp.load("Ylm_amplitude", emsbh_params.Ylm_amplitude, 0.);
        pp.load("Ylm_thickness", emsbh_params.Ylm_thickness, 3.);
        pp.load("Ylm_r0", emsbh_params.Ylm_r0 , 0.5 * L);


        // Do we want Weyl extraction, puncture tracking and constraint norm
        // calculation?
        pp.load("activate_extraction", activate_extraction, false);

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("num_mass_extraction_radii",
                mass_extraction_params.num_extraction_radii, 1);
        pp.load("mass_extraction_levels",
                mass_extraction_params.extraction_levels,
                mass_extraction_params.num_extraction_radii, 0);
        pp.load("mass_extraction_radii",
                mass_extraction_params.extraction_radii,
                mass_extraction_params.num_extraction_radii, 0.1);
        pp.load("num_points_phi_mass", mass_extraction_params.num_points_phi,
                2);
        pp.load("num_points_theta_mass",
                mass_extraction_params.num_points_theta, 4);
        pp.load("mass_extraction_center",
                mass_extraction_params.extraction_center,
                {0.0, 0.0});


        // Apparent Horizon stuff
        #ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, emsbh_params.bh_mass*0.5);
        #endif
        pp.load("AH_num_horizons", AH_num_horizons, 0);
        pp.load("AH_expect_merger", AH_expect_merger, 0);
        pp.load("horizon_centre_1", horizon_centre_1,
                {0.5 * L, 0.5 * L});
        pp.load("horizon_centre_2", horizon_centre_2,
                {0.5 * L, 0.5 * L});

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("num_mass_extraction_radii",
                mass_extraction_params.num_extraction_radii, 1);
        pp.load("mass_extraction_levels",
                mass_extraction_params.extraction_levels,
                mass_extraction_params.num_extraction_radii, 0);
        pp.load("mass_extraction_radii",
                mass_extraction_params.extraction_radii,
                mass_extraction_params.num_extraction_radii, 0.1);
        pp.load("num_points_phi_mass", mass_extraction_params.num_points_phi,
                2);
        pp.load("num_points_theta_mass",
                mass_extraction_params.num_points_theta, 4);
        pp.load("mass_extraction_center",
                mass_extraction_params.extraction_center,
                {0.5 * L, 0.5 * L});

        // Weyl extraction
        pp.load("activate_gw_extraction", activate_weyl_extraction, 0);

        // Em extraction stuff
        pp.load("activate_em_extraction", activate_pheyl_extraction, 0);

        // real scalar extraction stuff
        pp.load("activate_rs_extraction", activate_rs_extraction, 0);

        // mass charge scalar extraction stuff
        pp.load("activate_mq_extraction", activate_mq_extraction, 0);


        // Variables for outputting inf-norm
        pp.load("num_vars_inf_norm", num_vars_inf_norm, 0);
        pp.load("vars_inf_norm", vars_inf_norm, num_vars_inf_norm, 0);


        // Do we cant to calculate L2 norms of constraint violations
        pp.load("calculate_constraint_violations",
                calculate_constraint_violations, false);
    }

    void check_params()
    {}

    // tagging
    bool activate_extraction;

    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi, regrid_threshold_A;

    // extraction stuff
    int activate_weyl_extraction;
    int activate_pheyl_extraction;
    int activate_rs_extraction;
    int activate_mq_extraction;

    // // Layer regridding
    // bool m_do_layer_tagging;
    // double y_regrid_lim;

    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;

    // EMS BH stuff
    EMSBH_params_t emsbh_params;
    CouplingFunction::params_t coupling_function_params;

    double m_G_Newton;
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    // Vars for outputting inf-norms
    int num_vars_inf_norm;
    std::vector<int> vars_inf_norm;

#ifdef USE_AHFINDER
    double AH_initial_guess;
    int AH_num_horizons;
    int AH_expect_merger;
#endif
    std::array<double, CH_SPACEDIM> horizon_centre_1;
    std::array<double, CH_SPACEDIM> horizon_centre_2;

};
#endif /* SIMULATIONPARAMETERS_HPP_ */
