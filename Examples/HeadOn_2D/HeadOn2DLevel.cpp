/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "HeadOn2DLevel.hpp"
#include "BinaryPunctureTaggingCriterion.hpp"
#include "BoxLoops.hpp"
#include "CCZ4Cartoon.hpp"
#include "ComputePack.hpp"
#include "ConstraintsCartoon.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include "TraceARemovalCartoon.hpp"
#include "WeylExtraction.hpp"
#include "WeylOmScalar.hpp"

// Initial data
#include "HeadOn2D.hpp"
#include "TwoPuncturesInitialData.hpp"

// Reference connection
#include "GammaCartoonCalculator.hpp"

void HeadOn2DLevel::specificAdvance()
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

void HeadOn2DLevel::initialData()
{
    CH_TIME("HeadOn2DLevel::initialData");
    if (m_verbosity)
        pout() << "HeadOn2DLevel::initialData " << m_level << endl;

#ifdef USE_TWOPUNCTURES
    TwoPuncturesInitialData two_punctures_initial_data(
        m_dx, m_p.center, m_tp_amr.m_two_punctures);
    // Can't use simd with this initial data
    BoxLoops::loop(two_punctures_initial_data, //,TraceARemovalCartoon()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());
#else
    HeadOn2D binary(m_p.bh1_params, m_p.bh2_params, m_dx);
    // First set everything to zero (to avoid undefinded values)
    // then calculate initial data
    BoxLoops::loop(make_compute_pack(SetValue(0), binary), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS);
#endif

    // Needed only if Gamma's are non-zero (non-conformally flat ID)
    // BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new,
    // INCLUDE_GHOST_CELLS);
}

// Things to do before a plot level - need to calculate the Weyl scalars
void HeadOn2DLevel::prePlotLevel()
{
    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx), m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
}

void HeadOn2DLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                    const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side
    BoxLoops::loop(
        CCZ4Cartoon<>(m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
        a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

void HeadOn2DLevel::specificUpdateODE(GRLevelData &a_soln,
                                      const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void HeadOn2DLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                            const FArrayBox &current_state)
{
    std::array<double, 2> puncture_masses;
    std::vector<std::array<double, CH_SPACEDIM>> punctures;

    if (m_p.track_punctures)
    {
#ifdef USE_TWOPUNCTURES
        puncture_masses = {m_tp_amr.m_two_punctures.mm,
                           m_tp_amr.m_two_punctures.mp};
#else
        puncture_masses = {m_p.bh1_params.mass, m_p.bh2_params.mass};
#endif
        punctures = m_bh_amr.m_puncture_tracker.get_puncture_coords();
    }
    BoxLoops::loop(BinaryPunctureTaggingCriterion<FourthOrderDerivatives>(
                       m_dx, m_level, m_p.puncture_tag_max_levels,
                       m_p.extraction_params, punctures,
                       m_p.activate_extraction, m_p.track_punctures,
                       puncture_masses, m_p.bh_tagging_buffer,
                       m_p.puncture_tag_min_separation),
                   current_state, tagging_criterion); //, disable_simd());
}

void HeadOn2DLevel::specificPostTimeStep()
{
    CH_TIME("HeadOn2DLevel::specificPostTimeStep");
    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();

            BoxLoops::loop(WeylOmScalar(m_p.extraction_params.center, m_dx),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                m_gr_amr.m_interpolator->refresh();
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    // do puncture tracking on requested level
    if (m_p.track_punctures && m_level == m_p.puncture_tracking_level)
    {
        CH_TIME("PunctureTracking");
        // only do the write out for every coarsest level timestep
        int coarsest_level = 0;
        bool write_punctures = at_level_timestep_multiple(coarsest_level);
        m_bh_amr.m_puncture_tracker.execute_tracking(m_time, m_restart_time,
                                                     m_dt, write_punctures);
    }

#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        if (m_p.AH_set_origins_to_punctures && m_p.track_punctures)
        {
            m_bh_amr.m_ah_finder.set_origins(
                m_bh_amr.m_puncture_tracker.get_puncture_coords());
        }
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    }
#endif
}
