/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "BlackStringLevel.hpp"
#include "BlackStringTaggingCriterion.hpp"
#include "BoxLoops.hpp"
#include "CCZ4Cartoon.hpp"
#include "ComputePack.hpp"
#include "ConstraintsCartoon.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include "TraceARemovalCartoon.hpp"

// Initial data
#include "BlackString.hpp"

#include "GammaCartoonCalculator.hpp"

#ifdef USE_AHFINDER
#include "ApparentHorizon.hpp"
#endif

void BlackStringLevel::specificAdvance()
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

void BlackStringLevel::initialData()
{
    CH_TIME("BlackStringLevel::initialData");
    if (m_verbosity)
        pout() << "BlackStringLevel::initialData " << m_level << endl;

    BlackString blackstring(m_p.blackstring_params, m_dx, m_p.initial_lapse,
                            m_p.initial_shift); // Set up the compute class for
                                                // the BlackString initial data
    // First set everything to zero (to avoid undefinded values on constraints)
    // then calculate initial data
    BoxLoops::loop(make_compute_pack(SetValue(0), blackstring), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS);

    //    BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new,
    //    INCLUDE_GHOST_CELLS);
}

// Things to do before a plot level - need to calculate the Weyl scalars
void BlackStringLevel::prePlotLevel()
{
    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

void BlackStringLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side and set constraints to zero to avoid
    // undefined values
    BoxLoops::loop(make_compute_pack(SetValue(0, Interval(c_Ham, NUM_VARS - 1)),
                                     CCZ4Cartoon<>(m_p.ccz4_params, m_dx,
                                                   m_p.sigma, m_p.formulation)),
                   a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
} // new m_dt

void BlackStringLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void BlackStringLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(BlackStringTaggingCriterion(m_dx, m_p.chiTagMultiplier,
                                               m_p.KTagMultiplier,
                                               m_p.gammaTagMultiplier),
                   current_state, tagging_criterion);
}

void BlackStringLevel::specificPostTimeStep()
{
    CH_TIME("BlackStringLevel::specificPostTimeStep");
#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);

        if (m_bh_amr.m_ah_finder.get(0)->get_converged())
        {
            // warning if max level is not enough to cover the neck
            double minimum_y_ah = m_bh_amr.m_ah_finder.get(0)->get_min_F();
            double coarsest_dx = m_dx * pow(2., m_level);
            double required_dx = minimum_y_ah / m_p.AH_min_points;
            int required_level = ceil(log2(coarsest_dx / required_dx));

            pout() << "Horizon minimum at y = " << minimum_y_ah << std::endl;

            if (required_level > m_p.max_level)
            {
                pout() << "Add more levels. Current max_level = "
                       << m_p.max_level
                       << ". Recommended max_level = " << required_level
                       << " to reach horizon minimum at y = " << minimum_y_ah
                       << " with at least " << m_p.AH_min_points << " points."
                       << std::endl;
                if (m_p.AH_stop_if_problems)
                    MayDay::Error("Not enough levels. Problems! You asked me "
                                  "to stop. Check the pout's");
            }
        }
        else if (m_p.AH_stop_if_problems)
            MayDay::Error(
                "AH didn't converge. Problems! You asked me to stop.");
    }
#endif
}
