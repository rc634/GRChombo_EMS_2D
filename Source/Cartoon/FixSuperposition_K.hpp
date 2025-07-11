/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXSUPERPOSITION_K_HPP_
#define FIXSUPERPOSITION_K_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

class FixSuperposition_K
{
    // Only variables needed are metric
    template <class data_t> struct Vars
    {
        Tensor<2, data_t> h;
        Tensor<2, data_t> A;
        Tensor<1, data_t> shift;
        data_t hww;
        data_t Aww;
        data_t lapse;
        data_t chi;
        data_t K;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_h11, c_h22>(), h);

            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_A11, c_A22>(), A);

            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_shift1, c_shift2>(), shift);

            VarsTools::define_enum_mapping(mapping_function, c_hww, hww);

            VarsTools::define_enum_mapping(mapping_function, c_Aww, Aww);

            VarsTools::define_enum_mapping(mapping_function, c_lapse, lapse);

            VarsTools::define_enum_mapping(mapping_function, c_chi, chi);

            VarsTools::define_enum_mapping(mapping_function, c_K, K);
        }
    };

  protected:
    const FourthOrderDerivatives m_deriv;
    const double m_dx;
    std::array<double, CH_SPACEDIM> m_center;
    const double m_rapidity;
    const bool m_binary;

  public:
    FixSuperposition_K(
        double a_dx, std::array<double, CH_SPACEDIM> a_center, double a_rapidity, bool a_binary)
        : m_deriv(a_dx), m_dx(a_dx), m_center(a_center), m_rapidity(a_rapidity), m_binary(a_binary)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        //////////////////////////////
        // my stuff
        //////////////////////////////

        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        using namespace TensorAlgebra;

        // INPUT H IS ACTAULLY GAMMA HERE!!!!!!!
        // INPUTE A IS ACTUALLY (DT GAMMA / LAPSE) HERE!!!!
        // fix superposition error in physical metric and lapse
        double one_if_binary = 0.;
        if (m_binary)
        {
            one_if_binary = 1.;
        }

        // print statement that is only triggered for a couple central gridpoints
        if (coords.x*coords.x + coords.y*coords.y < 1.)
        {
          //std::cout << "one_if_binary = " << one_if_binary << std::endl;
        }

        double lapse = vars.lapse - one_if_binary;
        Tensor<2, data_t, 3> gij = {0.}; // spatial metric
        Tensor<2, data_t, 3> Aij = {0.};
        Tensor<2, data_t, 3> Kij = {0.};
        Tensor<2, data_t, 3> LieBeta_ij = {0.};
        Tensor<2, data_t, 3> dt_gamma_ij = {0.};

        // load physical metric from stored hij
        FOR2(i,j) gij[i][j] = vars.h[i][j];
        FOR1(i) gij[i][i] -= one_if_binary;
        gij[2][2] = vars.hww - one_if_binary;

        // // beta^2 for debug
        // double betasqr = 0.;
        // FOR2(i,j)
        // {
        //   betasqr += vars.shift[i] * vars.shift[j] * gij[i][j];
        // }


        // inverse physical metric
        const auto g_UU = compute_inverse_sym(gij);


        // load time derivative of physical metric from stored Aij
        FOR2(i,j) dt_gamma_ij[i][j] = 0.*vars.A[i][j];
        dt_gamma_ij[2][2] = 0.*vars.Aww;


        // d_i shift^j
        Tensor<2, data_t, 3> dshift_ij = {0.};
        // evaluate spatial partial derivatives of normal vector
        FOR2(i,j)
        {
           dshift_ij[i][j] = d1.shift[j][i];
        }
        dshift_ij[2][2] = vars.shift[1] / coords.y;


        // Eval Lie deriv
        // evaluate the first term of the lie derivative
        FOR3(i,j,k)
        {
            LieBeta_ij[i][j] += vars.shift[k] * d1.h[i][j][k];
        }
        FOR1(i) {
            LieBeta_ij[2][2] += vars.shift[i] * d1.hww[i];
        }
        // evaluate the other components of the lie derivative
        for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
        for (int k=0; k<3; k++){
            LieBeta_ij[i][j] += gij[k][j] * dshift_ij[i][k]
                              + gij[i][k] * dshift_ij[j][k];
        }}}


        // calcualte K_ij from the definition
        for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
          Kij[i][j] = - (0.5 / lapse) * (dt_gamma_ij[i][j] - LieBeta_ij[i][j]);
          // Kij[i][j] = dt_gamma_ij[i][j]; // debugging
          // Kij[i][j] = dshift_ij[i][j] - dt_gamma_ij[i][j];
        }}


        // take the trace of Kij -> K = gamma^{ij} K_{ij}
        double trK = 0.;
        for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
          trK += Kij[i][j] * g_UU[i][j];
          //trK += dt_gamma_ij[i][j] * g_UU[i][j]; //  this is for debugging
        }}


        // calculate conformal factor chi to decompose Kij into Aij
        double det_gamma = compute_determinant_sym(gij);
        double local_chi = pow(det_gamma,-1./3.);


        // calculate conformal traceless curvature Aij
        for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
          Aij[i][j] = local_chi * (Kij[i][j] - (1./3.) * trK * gij[i][j]);
          // Aij[i][j] = dt_gamma_ij[i][j] - d1.h[i][j][0]; // for debug
        }}


        // Store variables !
        // trace K
        // current_cell.store_vars(trK, c_K);
        // // Aij
        // current_cell.store_vars(Aij[2][2],c_Aww);
        // current_cell.store_vars(Aij[0][0],c_A11);
        // current_cell.store_vars(Aij[0][1],c_A12);
        // current_cell.store_vars(Aij[1][1],c_A22);

    }

};

#endif /* FIXSUPERPOSITION_K_HPP_ */
