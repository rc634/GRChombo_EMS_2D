/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(ISOTROPICBOOSTEDBBH_HPP_)
#error "This file should only be included through IsotropicBoostedBBH.hpp"
#endif

#ifndef ISOTROPICBOOSTEDBBH_IMPL_HPP_
#define ISOTROPICBOOSTEDBBH_IMPL_HPP_

// Compute the value of the initial vars on the grid
template <class data_t>
void IsotropicBoostedBBH::compute(Cell<data_t> current_cell) const
{
    Coordinates<data_t> coords1(current_cell, m_dx, m_bh1.m_bh_params.center);
    Coordinates<data_t> coords2(current_cell, m_dx, m_bh2.m_bh_params.center);

    Tensor<2, data_t, 3> gammaLL_1, gammaUU_1, KLL_1;
    data_t shift_x_1;

    Tensor<2, data_t, 3> gammaLL_2, gammaUU_2, KLL_2;
    data_t shift_x_2;

    // pre-computations
    m_bh1.computeMetric(gammaLL_1, gammaUU_1, shift_x_1, KLL_1, coords1);
    m_bh2.computeMetric(gammaLL_2, gammaUU_2, shift_x_2, KLL_2, coords2);

    // Superpose solutions
    Tensor<2, data_t, 3> gammaLL, gammaUU, KLL;
    data_t shift_x;

    data_t g_xx = gammaLL_1[0][0] + gammaLL_2[0][0] - 1.;
    data_t g_yy = gammaLL_1[1][1] + gammaLL_2[1][1] - 1.;
    data_t g_zz = gammaLL_1[2][2] + gammaLL_2[2][2] - 1.;

    gammaLL[0][0] = g_xx;
    gammaLL[1][1] = g_yy;
    gammaLL[2][2] = g_zz;
    gammaUU[0][0] = 1. / g_xx;
    gammaUU[1][1] = 1. / g_yy;
    gammaUU[2][2] = 1. / g_zz;

    shift_x = shift_x_1 + shift_x_2;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int m = 0; m < 3; m++)
            {
                for (int k = 0; k < 3; k++)
                {
                    KLL[i][j] += 0.5 * ((KLL_1[k][i] * gammaUU_1[k][m] +
                                         KLL_2[k][i] * gammaUU_2[k][m]) *
                                            gammaLL[j][m] +
                                        (KLL_1[k][j] * gammaUU_1[k][m] +
                                         KLL_2[k][j] * gammaUU_2[k][m]) *
                                            gammaLL[i][m]);
                }
            }
        }
    }

    IsotropicBoostedBH::Vars<data_t> vars;
    // compute vars
    IsotropicBoostedBH::computeVars(vars, gammaLL, gammaUU, shift_x, KLL);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* ISOTROPICBOOSTEDBBH_IMPL_HPP_ */
