// Last edited James Cook 3. August 2017
#ifndef OSCILLOTON_HPP_
#define OSCILLOTON_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
// #include "MatterCCZ4.hpp"
// #include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <array>
#include <vector>

// Class which creates Oscillotons (A selfgravitating scalar field) using
// external data-files. Information taken form gr-qc/0301105 - Numerical sutdies
// of Phi^2 Oscillotons
class Oscilloton
{
  public:
    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    struct params_t
    {
        // double amplitudeSF; // Amplitude of bump in initial SF bubble
        std::array<double, 3>
            centerSF; // Centre of perturbation in initial SF bubble
        std::array<double, 3>
            centerSF2;  // Centre of perturbation in initial SF bubble

        // double widthSF; // Width of bump in initial SF bubble
        // double r_zero;  // Position of bump relative to centre
        double vx;      // x-direction velocity of oscilloton
        double vx2;     // x-direction of 2nd Oscilloton
    };

    Oscilloton(params_t a_matter_params, double a_dx, double spacing);

  protected:
    double m_dx;
    double m_spacing;

    // Old Vectors for inputting from CSV
    vector<double> m_lapse_values;
    vector<double> m_grr_values;
    vector<double> m_Krr_values;
    vector<double> m_Phi_values;
    vector<double> m_Pi_values;

    // Calculated Values
    double m_omega; // Omega = sqrt(c0[final value])/a0[final value];
    double m_time;  // An initial 'unboosted' time
    double m_gamma; // Gamma Factor
    double m_gamma2;

    // lapse = (m_b - m_rs_alpha/rr)
    double m_rs_alpha; // r_shwarzchild calculated from alpha
    double m_b;        // r_shwarzchild modfied alpha

    // grr = (m_c - m_rs_g_rr/rr)^(-1)
    double m_rs_g_rr; // r_shwarzchild calculated from grr
    double m_c;       // r_shwarzchild modified grr

    // row and column max of input files
    int m_row_max; // Maximum rows in our arrays (Sorry Katy for using arrays
                   // and not tensors!)
    int m_column_max; // Maximum columns in our arrays (notice that phi has
                      // m_column_max-1)

    const static int m_file_length = 2150;

    // 2D arrays to hold the input data files
    double a202[m_file_length][8];        // a_j
    double c202[m_file_length][8];        // c_j
    double phi202[m_file_length][7];      // phi_j
    double d_a202_dr[m_file_length][8];   // dr a_j
    double d_c202_dr[m_file_length][8];   // dr c_j
    double d_phi202_dr[m_file_length][7]; // dr phi_j

    // To input the components of a, c and phi (and there respective
    // derivatives) into the linear Interpolator we need them to be a 1D array

    // a, c and phi components
    double a_rad[m_file_length], a0[m_file_length], a2[m_file_length],
        a4[m_file_length], a6[m_file_length], a8[m_file_length],
        a10[m_file_length], a12[m_file_length];
    double c_rad[m_file_length], c0[m_file_length], c2[m_file_length],
        c4[m_file_length], c6[m_file_length], c8[m_file_length],
        c10[m_file_length], c12[m_file_length];
    double phi_rad[m_file_length], phi1[m_file_length], phi3[m_file_length],
        phi5[m_file_length], phi7[m_file_length], phi9[m_file_length],
        phi11[m_file_length];

    // a, c and phi radial derivatives
    double d_a_rad_dr[m_file_length], d_a0_dr[m_file_length],
        d_a2_dr[m_file_length], d_a4_dr[m_file_length], d_a6_dr[m_file_length],
        d_a8_dr[m_file_length], d_a10_dr[m_file_length],
        d_a12_dr[m_file_length];
    double d_c_rad_dr[m_file_length], d_c0_dr[m_file_length],
        d_c2_dr[m_file_length], d_c4_dr[m_file_length], d_c6_dr[m_file_length],
        d_c8_dr[m_file_length], d_c10_dr[m_file_length],
        d_c12_dr[m_file_length];
    double d_phi_rad_dr[m_file_length], d_phi1_dr[m_file_length],
        d_phi3_dr[m_file_length], d_phi5_dr[m_file_length],
        d_phi7_dr[m_file_length], d_phi9_dr[m_file_length],
        d_phi11_dr[m_file_length];

    // New Linear Interpolator - 1D array with data (coords) to value
    //  template <class data_t>

    // New Linear Interpolator - 1D array with data (coords) to value
    // template <class data_t>
    //   double linear_interpolation_new2(Coordinates<data_t> coords, const
    //   double (*data_jt)) const ;

    // New Linear Interpolator - 1D array with data (coords) to value
    // template <class data_t>
    // double linear_interpolation_new3(Coordinates<data_t> coords, const double
    // (*data_jt)) const ;

  public:
    struct StarData
    {
        // Construct full 3D data here, cartoon down later
        double phi, Pi, lapse;
        double beta_U[3];
        double beta[3];
        Tensor<2, double, 3> h_cartesian;
        Tensor<2, double, 3> K_cartesian;
        Tensor<2, double, 3> h_cartesian_UU;

        // given params
        double StarX, StarY, StarZ;
        double vx;
    };

    const params_t m_params; //!< The matter params

    template <class data_t> void Spawn(StarData &DeStar) const;

    template <class data_t>
    double linear_interpolation_new(double vx, double t, double x, double y,
                                    double z, const double(*data_jt)) const;

    template <class data_t> void test2022(data_t type);
};

#include "Oscilloton.impl.hpp"

#endif /* OSCILLOTON_HPP_ */
