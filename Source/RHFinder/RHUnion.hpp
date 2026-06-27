#ifndef RHUNION_HPP_
#define RHUNION_HPP_

#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "RHSurf.hpp"
#include "UserVariables.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

// Owns all tracked horizons.
class RHUnion
{
  public:
    enum : int { NG = 2 }; // ghost cells on each side of each surface array
    int    m_num_stale_repeats = 5;       // stale repeats when err > thresh_high
    double m_courant           = 0.25;     // normal chase step size
    double m_thresh_high       = 0.0001;    // above: slow chase + stale repeats
    double m_thresh_close     = 0.0005;   // below: switch from fast chase to Newton
    double m_thresh_super_low  = 0.0000001; // below: surface converged, stop
    double m_newton_delta_f    = 1e-4;    // finite difference step for Jacobian assembly
    std::vector<RHSurf> m_surfaces;
    std::vector<std::ofstream> m_outfiles;
    std::vector<std::ofstream> m_ffiles;

    RHUnion() {}

    void setup(int a_num,
               const std::vector<double> &a_radii,
               const std::vector<double> &a_centres_x,
               const std::vector<int>    &a_num_points,
               const std::vector<int>    &a_levels,
               const std::vector<int>    &a_time_step_freqs,
               const std::vector<double> &a_newton_crits,
               const std::vector<double> &a_chase_speeds,
               const std::vector<double> &a_start_times,
               double                     a_restart_time = 0.)
    {
        m_surfaces.clear();
        m_outfiles.clear();
        m_outfiles.resize(a_num);
        m_ffiles.clear();
        m_ffiles.resize(a_num);
        for (int i = 0; i < a_num; ++i)
        {
            m_surfaces.emplace_back(a_num_points[i], NG, std::vector<double>{a_centres_x[i], 0.0}, i);
            m_surfaces.back().m_level          = a_levels[i];
            m_surfaces.back().m_time_step_freq = a_time_step_freqs[i];
            m_surfaces.back().m_newton_crit    = a_newton_crits[i];
            m_surfaces.back().m_chase_speed    = a_chase_speeds[i];
            m_surfaces.back().m_start_time     = a_start_times[i];
            std::fill(m_surfaces.back().m_f.begin(),
                      m_surfaces.back().m_f.end(), a_radii[i]);
            if (procID() == 0)
            {
                const bool is_restart = (a_restart_time > 0.);
                if (is_restart)
                    m_surfaces.back().prune(a_restart_time);
                const auto mode = is_restart
                    ? std::ios::out | std::ios::app
                    : std::ios::out | std::ios::trunc;

                m_outfiles[i].open("rh_surf_" + std::to_string(i) + ".dat", mode);
                if (!is_restart)
                    m_outfiles[i] << std::setw(16) << "#t"
                                  << std::setw(16) << "surf"
                                  << std::setw(16) << "x"
                                  << std::setw(16) << "<r>"
                                  << std::setw(16) << "Area"
                                  << std::setw(16) << "M_irr"
                                  << std::setw(16) << "P_x"
                                  << std::setw(16) << "<Theta+>"
                                  << std::setw(16) << "<Theta->"
                                  << std::setw(16) << "err"
                                  << std::setw(16) << "K0"
                                  << std::setw(16) << "K2"
                                  << std::setw(16) << "K4"
                                  << std::setw(16) << "KP2"
                                  << std::setw(16) << "KP4"
                                  << std::setw(16) << "eq_perim"
                                  << std::setw(16) << "pol_perim"
                                  << std::setw(16) << "Q_charge"
                                  << std::setw(16) << "M_total"
                                  << std::setw(8)  << "mode"
                                  << std::endl;

                m_ffiles[i].open("rh_f" + std::to_string(i) + ".dat", mode);
                if (!is_restart)
                {
                    m_ffiles[i] << std::setw(16) << "#t";
                    for (int j = 0; j < a_num_points[i]; ++j)
                        m_ffiles[i] << std::setw(16) << ("f[" + std::to_string(j) + "]");
                    m_ffiles[i] << std::endl;
                }
            }
        }

        if (procID() == 0)
        {
            pout() << "\nRHUnion::Construction  Initial surfaces are:" << std::endl;
            for (const auto &surf : m_surfaces)
                surf.hello();
        }
    }

    // Set EMS coupling parameters on all surfaces so Q_charge() uses the right coupling.
    void set_coupling_params(double alpha, double f0, double f1, double f2)
    {
        for (auto &surf : m_surfaces)
        {
            surf.m_alpha = alpha;
            surf.m_f0    = f0;
            surf.m_f1    = f1;
            surf.m_f2    = f2;
        }
    }

    void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }

    // Interpolate all needed fields at interior surface points in one combined
    // query, then store results in each surface's arrays and fill ghost cells.
    void interpolate_fields()
    {
        int total_pts = 0;
        for (const auto &surf : m_surfaces) {
            total_pts += surf.m_n;
        }

        std::vector<double> qx(total_pts), qy(total_pts);
        int offset = 0;
        for (const auto &surf : m_surfaces)
        {
            for (int i = 0; i < surf.m_n; ++i)
            {
                int ii = i + surf.m_NG; // non-ghost array index
                qx[offset + i] = surf.m_centre[0] + surf.m_f[ii] * std::cos(surf.m_theta[ii]);
                qy[offset + i] = surf.m_centre[1] + surf.m_f[ii] * std::sin(surf.m_theta[ii]);
            }
            offset += surf.m_n;
        }

        // output buffers — values
        std::vector<double> chi(total_pts), K(total_pts);
        std::vector<double> A11(total_pts), A12(total_pts), A22(total_pts), Aww(total_pts);
        std::vector<double> h11(total_pts), h12(total_pts), h22(total_pts), hww(total_pts);
        // output buffers — derivatives of chi and conformal metric
        std::vector<double> dx_chi(total_pts), dy_chi(total_pts);
        std::vector<double> dx_h11(total_pts), dy_h11(total_pts);
        std::vector<double> dx_h12(total_pts), dy_h12(total_pts);
        std::vector<double> dx_h22(total_pts), dy_h22(total_pts);
        std::vector<double> dx_hww(total_pts), dy_hww(total_pts);
        // EM and scalar field buffers
        std::vector<double> Ex(total_pts), Ey(total_pts), phi(total_pts);

        InterpolationQuery query(total_pts);
        query.setCoords(0, qx.data()).setCoords(1, qy.data());

        query.addComp(c_chi, chi.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_K,   K.data(),   Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_A11, A11.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_A12, A12.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_A22, A22.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_Aww, Aww.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_h11, h11.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_h12, h12.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_h22, h22.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_hww, hww.data(), Derivative::LOCAL, VariableType::evolution);

        query.addComp(c_chi, dx_chi.data(), Derivative::dx, VariableType::evolution);
        query.addComp(c_chi, dy_chi.data(), Derivative::dy, VariableType::evolution);
        query.addComp(c_h11, dx_h11.data(), Derivative::dx, VariableType::evolution);
        query.addComp(c_h11, dy_h11.data(), Derivative::dy, VariableType::evolution);
        query.addComp(c_h12, dx_h12.data(), Derivative::dx, VariableType::evolution);
        query.addComp(c_h12, dy_h12.data(), Derivative::dy, VariableType::evolution);
        query.addComp(c_h22, dx_h22.data(), Derivative::dx, VariableType::evolution);
        query.addComp(c_h22, dy_h22.data(), Derivative::dy, VariableType::evolution);
        query.addComp(c_hww, dx_hww.data(), Derivative::dx, VariableType::evolution);
        query.addComp(c_hww, dy_hww.data(), Derivative::dy, VariableType::evolution);
        query.addComp(c_Ex,  Ex.data(),     Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_Ey,  Ey.data(),     Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_phi, phi.data(),    Derivative::LOCAL, VariableType::evolution);

        m_interpolator->interp(query);

        offset = 0;
        for (auto &surf : m_surfaces)
        {
            for (int i = 0; i < surf.m_n; ++i)
            {
                int ii = surf.m_NG + i; // non-ghost array index
                surf.m_chi[ii]    = chi[offset + i];
                surf.m_K[ii]      = K[offset + i];
                surf.m_A11[ii]    = A11[offset + i];
                surf.m_A12[ii]    = A12[offset + i];
                surf.m_A22[ii]    = A22[offset + i];
                surf.m_Aww[ii]    = Aww[offset + i];
                surf.m_h11[ii]    = h11[offset + i];
                surf.m_h12[ii]    = h12[offset + i];
                surf.m_h22[ii]    = h22[offset + i];
                surf.m_hww[ii]    = hww[offset + i];
                surf.m_dx_chi[ii] = dx_chi[offset + i];
                surf.m_dy_chi[ii] = dy_chi[offset + i];
                surf.m_dx_h11[ii] = dx_h11[offset + i];
                surf.m_dy_h11[ii] = dy_h11[offset + i];
                surf.m_dx_h12[ii] = dx_h12[offset + i];
                surf.m_dy_h12[ii] = dy_h12[offset + i];
                surf.m_dx_h22[ii] = dx_h22[offset + i];
                surf.m_dy_h22[ii] = dy_h22[offset + i];
                surf.m_dx_hww[ii] = dx_hww[offset + i];
                surf.m_dy_hww[ii] = dy_hww[offset + i];
                surf.m_Ex[ii]     = Ex[offset + i];
                surf.m_Ey[ii]     = Ey[offset + i];
                surf.m_phi[ii]    = phi[offset + i];
            }
            surf.fill_all_ghosts();
            offset += surf.m_n;
        }
    }

    // One finder step for all surfaces assigned to a_level.
    // Flow: dormancy check → interpolate → classify → chase loop → output.
    void update(double a_time, int a_level)
    {
        // skip entirely if no surface runs on this AMR level
        bool any = false;
        for (const auto &s : m_surfaces)
            if (s.m_level == a_level) { any = true; break; }
        if (!any) return;

        const auto t_start = std::chrono::steady_clock::now();

        // reset any surface whose mean radius has drifted out of bounds
        for (auto &surf : m_surfaces)
        {
            if (surf.m_dead || surf.m_state == RHSurf::SolverState::DORMANT) continue;
            const double r = surf.average_f();
            if (r <= surf.m_rmin || r >= surf.m_rmax)
            {
                pout() << "RHFinder: surface " << surf.m_index
                       << " out of bounds (mean_r=" << r << "), resetting to r=5\n";
                surf.reset(5.0);
            }
        }

        // update dormancy: surfaces before their start time stay/become dormant;
        // surfaces that have just passed their start time wake to FAR for classification
        for (auto &surf : m_surfaces)
        {
            if (surf.m_dead) continue;
            if (a_time < surf.m_start_time)
                surf.m_state = RHSurf::SolverState::DORMANT;
            else if (surf.m_state == RHSurf::SolverState::DORMANT)
                surf.m_state = RHSurf::SolverState::FAR;
        }

        // pull fresh grid data for all surfaces in one MPI round-trip
        m_interpolator->refresh();
        interpolate_fields();

        // re-centre all surfaces using the fresh field data
        for (auto &surf : m_surfaces)
            surf.re_centre();

        const int n = (int)m_surfaces.size();

        // measure pre-chase expansion error and assign solver regime;
        // this regime is held fixed for the entire chase below
        std::vector<double> errs(n);
        for (int k = 0; k < n; ++k)
        {
            auto &surf = m_surfaces[k];
            if (surf.m_state == RHSurf::SolverState::DORMANT) { errs[k] = 0.; continue; }
            errs[k] = surf.expansion_error();
            if      (errs[k] <= m_thresh_super_low) surf.m_state = RHSurf::SolverState::FOUND;
            else if (errs[k] <= m_thresh_high)      surf.m_state = RHSurf::SolverState::CLOSE;
            else                                    surf.m_state = RHSurf::SolverState::FAR;
        }

        // precompute per-surface courants; max_iters is the longest active surface's quota
        std::vector<double> courants(n);
        int max_iters = 0;
        for (int k = 0; k < n; ++k)
        {
            courants[k] = m_courant * m_surfaces[k].m_chase_speed;
            if (!m_surfaces[k].m_dead &&
                m_surfaces[k].m_state != RHSurf::SolverState::DORMANT &&
                m_surfaces[k].m_level == a_level &&
                errs[k] > m_thresh_super_low)
                max_iters = std::max(max_iters, m_surfaces[k].m_time_step_freq);
        }

        // chase — step-first so all active surfaces advance in lockstep and share
        // the same interpolated field snapshot at every iteration.
        // FAR surfaces take extra stale steps to close the gap faster.
        for (int i = 0; i < max_iters; ++i)
        {
            for (int k = 0; k < n; ++k)
            {
                auto &surf = m_surfaces[k];
                if (surf.m_level != a_level || surf.m_dead) continue;
                if (surf.m_state == RHSurf::SolverState::DORMANT) continue;
                if (errs[k] <= m_thresh_super_low) continue; // FOUND: nothing to do
                if (i >= surf.m_time_step_freq) continue;    // surface finished its quota

                try
                {
                    surf.chase_step(courants[k]);
                    if (errs[k] > m_thresh_high) // FAR: stale repeats to close gap faster
                        for (int j = 0; j < m_num_stale_repeats; ++j)
                            surf.chase_step(courants[k]);
                }
                catch (const std::exception &e)
                {
                    surf.m_dead = true;
                    pout() << "\nRHFinder: surface " << k << " marked dead: "
                           << e.what() << "\n";
                }
            }
            interpolate_fields(); // one batched MPI call refreshes all surfaces

            // early exit if every active surface on this level is converged
            bool all_found = true;
            for (int k = 0; k < n; ++k)
            {
                const auto &surf = m_surfaces[k];
                if (surf.m_level != a_level || surf.m_dead) continue;
                if (surf.m_state == RHSurf::SolverState::DORMANT) continue;
                if (surf.expansion_error() > m_thresh_super_low) { all_found = false; break; }
            }
            if (all_found) break;
        }

        // Newton polish — activated per surface once err drops below m_newton_crit.
        // Reverts if Newton diverges.
        // for (int k = 0; k < n; ++k)
        // {
        //     auto &surf = m_surfaces[k];
        //     if (surf.m_level != a_level || surf.m_dead) continue;
        //     if (errs[k] <= m_thresh_super_low) continue;
        //     if (surf.m_newton_crit <= 0.0) continue;

        //     try
        //     {
        //         const double err_post = surf.expansion_error();
        //         if (err_post < surf.m_newton_crit)
        //         {
        //             const std::vector<double> f_before_newton = surf.m_f;

        //             for (int ns = 0; ns < 10; ++ns)
        //             {
        //                 interpolate_fields();
        //                 surf.banded_newton_step(m_newton_delta_f, m_courant);
        //             }
        //             interpolate_fields();

        //             const double err_after_newton = surf.expansion_error();
        //             if (err_after_newton > err_post)
        //             {
        //                 surf.m_f = f_before_newton;
        //                 surf.fill_all_ghosts();
        //                 interpolate_fields();
        //                 pout() << "RHFinder: Newton worsened surface " << k
        //                        << " (" << err_post << " -> " << err_after_newton
        //                        << "), reverting.\n";
        //             }
        //         }
        //     }
        //     catch (const std::exception &e)
        //     {
        //         surf.m_dead = true;
        //         pout() << "\nRHFinder: surface " << k << " marked dead: "
        //                << e.what() << "\n";
        //     }
        // }

        // re-evaluate states post-chase so output and display reflect the final error
        for (int k = 0; k < n; ++k)
        {
            auto &surf = m_surfaces[k];
            if (surf.m_dead) continue;
            if (surf.m_state == RHSurf::SolverState::DORMANT) continue;
            const double err = surf.expansion_error();
            if      (err <= m_thresh_super_low) surf.m_state = RHSurf::SolverState::FOUND;
            else if (err <= m_thresh_high)      surf.m_state = RHSurf::SolverState::CLOSE;
            else                                surf.m_state = RHSurf::SolverState::FAR;
        }

        // write diagnostics (rhout) and full f(theta) profile (fout) for all surfaces
        rhout(a_time);
        fout(a_time);

        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_start).count();

        if (procID() == 0)
        {
            pout() << "\nRHFinder::update  t = " << std::fixed << std::setprecision(4) << a_time
                   << "  (" << std::setprecision(3) << elapsed << " s)";
            for (const auto &surf : m_surfaces)
                surf.fancy_hello();
        }
    }

    void rhout(double a_time)
    {
        if (procID() != 0) return;
        for (int i = 0; i < (int)m_surfaces.size(); ++i)
        {
            if (m_surfaces[i].m_state == RHSurf::SolverState::DORMANT) continue;
            m_surfaces[i].print_diagnostics(m_outfiles[i], a_time);
        }
    }

    void fout(double a_time)
    {
        if (procID() != 0) return;
        for (int i = 0; i < (int)m_surfaces.size(); ++i)
        {
            if (m_surfaces[i].m_state == RHSurf::SolverState::DORMANT) continue;
            m_surfaces[i].print_f(m_ffiles[i], a_time);
        }
    }

    // Remove rows with t > a_restart_time from all output files.
    // Call this after setup() (so surface indices exist), then re-call setup()
    // with a_append = true so subsequent writes continue from the pruned tail.
    void prune_all(double a_restart_time)
    {
        if (procID() != 0) return;
        for (const auto &surf : m_surfaces)
            surf.prune(a_restart_time);
    }

    void hello() const
    {
        for (const auto &surf : m_surfaces)
        {
            surf.hello();
        }
    }

  private:
    AMRInterpolator<Lagrange<4>> *m_interpolator = nullptr;
};

#endif /* RHUNION_HPP_ */
