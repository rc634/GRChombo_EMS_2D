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
#include "BlackString.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // Initial data
        pp.load("pert_freq", blackstring_params.pert_freq);
        pp.load("pert_amp", blackstring_params.pert_amp);
        pp.load("regFrac", blackstring_params.regFrac);
        pp.load("initial_lapse", initial_lapse);
        pp.load("initial_shift", initial_shift);
        pp.load("center", blackstring_params.center);
        pp.load("Lhor", blackstring_params.Lhor);
        pp.load("N1", blackstring_params.N1);
        pp.load("N2", blackstring_params.N2);
        pp.load("L", blackstring_params.L);

        // Tagging
        pp.load("chiTagMultiplier", chiTagMultiplier);
        pp.load("KTagMultiplier", KTagMultiplier);
        pp.load("gammaTagMultiplier", gammaTagMultiplier);

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess);
        pp.load("AH_min_points", AH_min_points);
        pp.load("AH_stop_if_problems", AH_stop_if_problems);
#endif
    }

    // tagging
    double chiTagMultiplier, KTagMultiplier, gammaTagMultiplier;
    // Initial data
    int initial_lapse, initial_shift;

    // Collection of parameters necessary for the CCZ4 RHS
    BlackString::params_t blackstring_params;

#ifdef USE_AHFINDER
    double AH_initial_guess;
    int AH_min_points;
    bool AH_stop_if_problems;
#endif
};
#endif /* SIMULATIONPARAMETERS_HPP_ */
