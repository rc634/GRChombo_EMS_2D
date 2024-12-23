/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __    _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to LICENSE, in Chombo's root directory.
 */
#endif

// General includes:
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>

#include "parstream.H" //Gives us pout()
using std::endl;
#include "BHAMR.hpp"

#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "AMRInterpolator.hpp"
#include "ApparentHorizonTest3DLevel.hpp"
#include "Lagrange.hpp"
#include "SmallDataIO.hpp"

#ifdef USE_AHFINDER
#include "AHFinder.hpp"
#endif

int runApparentHorizonTest3D(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    GRParmParse pp(0, argv + argc, NULL, in_string.c_str());
    SimulationParameters sim_params(pp);

    BHAMR bh_amr;
    DefaultLevelFactory<ApparentHorizonTest3DLevel> ah_test_level_fact(
        bh_amr, sim_params);
    setupAMRObject(bh_amr, ah_test_level_fact);

    int status = 0;

#ifdef USE_AHFINDER
    // Set up interpolator and PETSc subcommunicator when AH extraction is
    // active
    AMRInterpolator<Lagrange<4>> interpolator(
        bh_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    bh_amr.set_interpolator(&interpolator);

    AHSurfaceGeometry sph(sim_params.center);
    bh_amr.m_ah_finder.add_ah(sph, sim_params.initial_guess,
                              sim_params.AH_params);

    if (!bh_amr.m_ah_finder.get(0)->get_converged())
        status = 3;
    else
    {

        // get area from file to determine status
        auto stats = SmallDataIO::read("stats_AH1.dat");
        if (stats.size() == 0)
            status = 2;
        else
        {
            double mass = stats[3][0];
            double spin = stats[5][0];
            double spin_x = stats[6][0] * mass;

            double error_perc =
                fabs(1. - mass / sim_params.kerr_params.mass) * 100;
            if (sim_params.kerr_params.spin != 0.)
            {
                error_perc +=
                    fabs(1. - spin / sim_params.kerr_params.spin) * 100;
                error_perc +=
                    fabs(1. - spin_x / sim_params.kerr_params.spin) * 100;
            }
            pout() << "error = " << error_perc << "%" << std::endl;
            status |= (error_perc > 0.1) ||
                      +std::isnan(
                          error_perc); // accept 0.1% error in area calculation
        }
    }
#endif

    return status;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runApparentHorizonTest3D(argc, argv);

    if (status == 0)
    {
        std::cout << "ApparentHorizon3D test passed." << endl;
        pout() << "ApparentHorizon3D test passed." << endl;
    }
    else
    {
        std::cout << "ApparentHorizon3D test FAILED with return code " << status
                  << endl;
        pout() << "ApparentHorizon3D test FAILED with return code " << status
               << endl;
    }

    mainFinalize();
    return status;
}
