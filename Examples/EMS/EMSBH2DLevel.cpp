/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
#include <fstream>
#include <mutex>
#include <memory>

#include "EMSBH2DLevel.hpp"

// General headers
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"

// For tagging cells
#include "EMSExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "CCZ4Cartoon.hpp"
#include "MovingPunctureGauge.hpp"
#include "ConstraintsCartoon.hpp"
#include "SetValue.hpp"
#include "TraceARemovalCartoon.hpp"

// EMS includes
#include "EMSBH_read.hpp"
#include "EMSCouplingFunction.hpp"

// For GW extraction
#include "WeylExtraction.hpp"
#include "WeylOmScalar.hpp"

// For MQ extraction
#include "CrudeMassChargeExtraction.hpp"

// FOR EM extraction_level
#include "Pheyl2.hpp"
#include "PheylExtraction.hpp"

// For Scalar Field integrator
#include "RealScalarExtraction.hpp"

// For analysis outputs
#include "EMSCartoonLorentzScalars.hpp"
#include "SmallDataIO.hpp"

// For EM Tensor in 3D
// #include "EMTensor.hpp"

// other
#include "ADMQuantities.hpp"
#include "ADMQuantitiesExtraction.hpp"
#include "GammaCartoonCalculator.hpp"

void EMSBH2DLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance: "),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void EMS2DBHLevel::initialData()
{
    CH_TIME("EMSBH2DLevel::initialData");
    if (m_verbosity)
        pout() << "EMSBH2DLevel::initialData " << m_level << endl;

    // First initalise a EMSBH object
    EMSBH_read emdbh(m_p.emdbh_params, m_p.coupling_function_params,
                         m_p.G_Newton, m_dx, m_verbosity);

   if (m_verbosity)
       pout() << "EMSBH2DLevel::initialData - Compute_1d_solution " << m_level << endl;
    // the max radius the code might need to calculate out to is L*sqrt(3)
    emdbh.compute_1d_solution(2. * m_p.L);


    if (m_verbosity)
        pout() << "EMSBH2DLevel::initialData - Interpolate to 3D grid " << m_level << endl;
    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then  initial conditions for EMSBH
    BoxLoops::loop(make_compute_pack(SetValue(0.0), emdbh), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS, disable_simd());


    if (m_verbosity)
        pout() << "EMSBH2DLevel::initialData - GammaCalc " << m_level << endl;
    BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts();
}

// Things to do before a plot level - need to calculate the Weyl scalars
void EMSBH2DLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    // ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(
        make_compute_pack(

            WeylOmScalar(m_p.extraction_params.center, m_dx),

            Constraints<CouplingFunction>(m_dx, potential, m_p.m_G_Newton)

            //EMSCartoonLorentzScalars<EinsteinMaxwellScalarFieldWithCoupling>(m_dx,
                                   //m_p.extraction_params.extraction_center,
                                             //m_p.coupling_function_params),

            //EMTensor<EinsteinMaxwellScalarFieldWithCoupling>(
                //emd_field, m_dx, c_rho,
                //Interval(c_s1, c_s3), Interval(c_s11, c_s33)),

            //Pheyl2(m_p.extraction_params.extraction_center,
                              //m_p.coupling_function_params, m_dx)
                                                                 ),

      m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void EMSBH2DLevel::specificEvalRHS(GRLevelData &a_soln,
                                         GRLevelData &a_rhs,
                                         const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side
    CouplingFunction coupling_function(m_p.coupling_function_params);
    CCZ4Cartoon<MovingPunctureGauge, FourthOrderDerivatives, CouplingFunction>(
        m_p.ccz4_params, m_dx, m_p.sigma, potential, m_p.m_G_Newton,
        m_p.formulation);
    SetValue set_analysis_vars_zero(0.0, Interval(c_Xi + 1, NUM_VARS - 1));
    auto compute_pack =
        make_compute_pack(my_ccz4_cartoon, set_analysis_vars_zero);
    BoxLoops::loop(compute_pack, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

void EMSBH2DLevel::specificUpdateODE(GRLevelData &a_soln,
                                           const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}



// probably delete this ??
// void EMSBH2DLevel::preTagCells(double &a_maximum)
// {
//     // Find at which radius the maximum rho value is so we can use it in the
//     // tagging criterion class
//     // Don't do this when we're not restarting from checkpoint at t=0 though. It
//     // crashes the code because it starts requesting data from higher levels
//     // before they are created.
//     if (!m_p.restart_from_checkpoint && m_time == 0.0)
//     {
//         a_maximum = 0.0;
//     }
//     else
//     {
//         AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
//         RealVect max_rho_index = amr_reductions.maxIndex(c_rho);
//         double r2 = .0000001 + max_rho_index[0] * max_rho_index[0] +
//                     max_rho_index[1] * max_rho_index[1];
// #if CH_SPACEDIM == 3
//         r2 += max_rho_index[2] * max_rho_index[2];
// #endif
//         a_maximum = sqrt(r2);
//     }
// }
//



void EMSBH2DLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                                 const FArrayBox &current_state)
{
    BoxLoops::loop(EMDExtractionTaggingCriterion(
                     m_dx, m_level, m_p.mass_extraction_params,
                     m_p.regrid_threshold_A, m_p.regrid_threshold_phi,
                     m_p.regrid_threshold_chi), current_state,
                                                tagging_criterion);

}

void EMSBH2DLevel::specificPostTimeStep()
{
    CH_TIME("EMSBH2DLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main

    fillAllGhosts();

    Potential potential(m_p.potential_params);
    // ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(WeylOmScalar(m_p.extraction_params.center, m_dx),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);
    BoxLoops::loop(Constraints<Potential>(m_dx, potential, m_p.m_G_Newton),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);


#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        // Hack: if avg radius of found AH is negative, we reset initial guess
        double current_avg_AH_radius =  m_bh_amr.m_ah_finder.get(0)->get_ave_F();
        if (current_avg_AH_radius < 0.)
        {
            m_bh_amr.m_ah_finder.get(0)->solver.reset_initial_guess();
            pout() << "AHFinder: resetting initial guess as avg radius is "
                      "negative."
                   << endl;
        }
        else if (current_avg_AH_radius > 10.)
        {
            m_bh_amr.m_ah_finder.get(0)->solver.reset_initial_guess();
            pout() << "AHFinder: resetting initial guess as avg radius is "
                      "too large."
                   << endl;
        }
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    }
#endif

    if (m_p.activate_extraction == 1 &&
       at_level_timestep_multiple(m_p.extraction_params.min_extraction_level()))
    {
        // Do the extraction on the min extraction level
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            if (m_verbosity)
            {
                pout() << "EMSBH2DLevel::specificPostTimeStep:"
                          " Extracting gravitational waves." << endl;
            }


            // Refresh the interpolator and do the interpolation
            m_gr_amr.m_interpolator->refresh();
            WeylExtraction gw_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gw_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }

    // if (at_level_timestep_multiple(0))
    // {
    //     BoxLoops::loop(NoetherCharge(), m_state_new, m_state_diagnostics,
    //               EXCLUDE_GHOST_CELLS);
    // }

    if (m_level == 0)
    {
        bool first_step = (m_time == 0.);
        AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
        double L2_Ham = amr_reductions.norm(c_Ham, 2, true);
        double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom2), 2, true);
        SmallDataIO constraints_file(m_p.data_path + "constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
        {
            constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
        }
        constraints_file.write_time_data_line({L2_Ham, L2_Mom});

        int adm_min_level = 0;
        bool calculate_adm = at_level_timestep_multiple(adm_min_level);
        if (calculate_adm)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
            double M_ADM = amr_reductions.sum(c_rho_ADM);
            SmallDataIO M_ADM_file(m_p.data_path + "M_ADM", m_dt, m_time,
                                m_restart_time, SmallDataIO::APPEND,
                                first_step);
            M_ADM_file.remove_duplicate_time_data();
            if (first_step)
            {
                M_ADM_file.write_header_line({"M_ADM"});
            }
            M_ADM_file.write_time_data_line({M_ADM});
        }

        // double noether_charge = amr_reductions.sum(c_N);
        // SmallDataIO noether_charge_file("NoetherCharge", m_dt, m_time,
        //                                 m_restart_time,
        //                                 SmallDataIO::APPEND,
        //                                 first_step);
        // noether_charge_file.remove_duplicate_time_data();
        // if (m_time == 0.)
        // {
        //     noether_charge_file.write_header_line({"Noether Charge"});
        // }
        // noether_charge_file.write_time_data_line({noether_charge});

        double min_chi = amr_reductions.min(c_chi);
        SmallDataIO min_chi_file(m_p.data_path + "min_chi",
                                        m_dt, m_time, m_restart_time,
                                        SmallDataIO::APPEND, first_step);
        min_chi_file.remove_duplicate_time_data();
        if (first_step)
        {
            min_chi_file.write_header_line({"min_chi"});
        }
        min_chi_file.write_time_data_line({min_chi});


        // Do Electromagnetic Radiation Integration
        if (m_p.activate_pheyl_extraction == 1 &&
            at_level_timestep_multiple(
                m_p.extraction_params.min_extraction_level()))
        {
            CH_TIME("EMDBS2DLevel::doAnalysis::EMRadExtraction");

            fillAllGhosts();
            // CouplingFunction coupling_function(m_p.coupling_function_params);
            // EinsteinMaxwellScalarFieldWithCoupling emd_field(coupling_function);

            // fill grid with pheyl2 im and real components
            auto pheyl2_compute_pack = make_compute_pack(
                Pheyl2(m_p.extraction_params.extraction_center,
                                            m_p.coupling_function_params, m_dx));
            BoxLoops::loop(pheyl2_compute_pack, m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);


            // Do the extraction on the min extraction level
            if (m_level == m_p.extraction_params.min_extraction_level())
            {
                if (m_verbosity)
                {
                    pout() << "EMSBH2DLevel::specificPostTimeStep:"
                              " Extracting electromagnetic waves."
                           << endl;
                }

                // Refresh the interpolator and do the interpolation
                m_bh_amr.m_interpolator->refresh();
                PheylExtraction em_extraction(m_p.pheyl2_extraction_params,
                                             m_dt, m_time,
                                             first_step, m_restart_time);
                em_extraction.execute_query(m_bh_amr.m_interpolator);
            }
    }
}