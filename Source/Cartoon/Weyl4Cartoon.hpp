/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 // protective wrapper to stop double definition
 #ifndef TETRAD_ALREADY_DEFINED_
 #define TETRAD_ALREADY_DEFINED_
 //! Struct for the null tetrad
 template <class data_t> struct Tetrad_t
 {
     Tensor<1, data_t, 3> u; //!< the vector u^i
     Tensor<1, data_t, 3> v; //!< the vector v^i
     Tensor<1, data_t, 3> w; //!< the vector w^i
 };
 #endif



#ifndef WEYL4_HPP_
#define WEYL4_HPP_

#include "CCZ4CartoonVars.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "simd.hpp"
#include <array>

//! Struct for the E and B fields
template <class data_t> struct EBFields_t
{
    Tensor<2, data_t, 3> E; //!< Electric component of Weyltensor
    Tensor<2, data_t, 3> B; //!< Magnetic component of Weyltensor
};



//! Struct for the Newman Penrose scalar
template <class data_t> struct NPScalar_t
{
    data_t Real; // Real component
    data_t Im;   // Imaginary component
};

template <class data_t> struct CCZ4vars3D
{
    data_t lapse;
    data_t K;
    data_t chi;
    Tensor<1, data_t, 3> B;
    Tensor<1, data_t, 3> shift;
    Tensor<1, data_t, 3> Gamma;
    Tensor<2, data_t, 3> h;
    Tensor<2, data_t, 3> A;

    Tensor<1, data_t, 3> d1_lapse;
    Tensor<1, data_t, 3> d1_K;
    Tensor<1, data_t, 3> d1_chi;
    Tensor<2, data_t, 3> d1_B;
    Tensor<2, data_t, 3> d1_shift;
    Tensor<2, data_t, 3> d1_Gamma;
    Tensor<3, data_t, 3> d1_h;
    Tensor<3, data_t, 3> d1_A;

    Tensor<2, data_t, 3> d2_chi;
    //    Tensor<2,Tensor<2,data_t,3>,3> d2_h;
    std::array<std::array<std::array<std::array<data_t, 4>, 4>, 4>, 4> d2_h;
};

//!  Calculates the Weyl4 scalar for spacetimes without matter content
/*!
   This class calculates the Weyl4 scalar real and im parts using definitions
   from Alcubierres book "Introduction to 3+1 Numerical Relativity". We use a
   decomposition of the Weyl tensor in electric and magnetic parts, which then
   is used together with the tetrads defined in "gr-qc/0104063" to calculate the
   Weyl4 scalar.
*/
class Weyl4
{
  public:
    // Use the variable definitions containing the needed quantities
    template <class data_t> using Vars = CCZ4CartoonVars::VarsWithGauge<data_t>;
    template <class data_t>
    using Diff2Vars = CCZ4CartoonVars::Diff2VarsNoGauge<data_t>;

    //! Constructor of class Weyl4
    /*!
        Takes in the centre for the calculation of the tetrads, and grid spacing
    */
    Weyl4(const std::array<double, CH_SPACEDIM> a_center, const double a_dx)
        : m_center(a_center), m_dx(a_dx), m_deriv(a_dx)
    {
    }

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
    const double m_dx;                              //!< the grid spacing
    const FourthOrderDerivatives m_deriv; //!< for calculating derivs of vars

    //! Calculation of Weyl_4 scalar
    template <class data_t>
    NPScalar_t<data_t> compute_Weyl4(const EBFields_t<data_t> &ebfields,
                                     const CCZ4vars3D<data_t> &vars,
                                     const Coordinates<data_t> &coords) const;

    //! Calculation of the tetrads
    template <class data_t>
    Tetrad_t<data_t>
    compute_null_tetrad(const CCZ4vars3D<data_t> &vars,
                        const Coordinates<data_t> &coords) const;

    //! Calulation of the decomposition of the Weyl tensor in Electric and
    //! Magnetic fields
    template <class data_t>
    EBFields_t<data_t>
    compute_EB_fields(const CCZ4vars3D<data_t> &vars,
                      const Coordinates<data_t> &coords) const;

    template <class data_t>
    CCZ4vars3D<data_t> load_ccz4(const Vars<data_t> &vars,
                                 const Vars<Tensor<1, data_t>> &d1,
                                 const Diff2Vars<Tensor<2, data_t>> &d2,
                                 const Coordinates<data_t> &coords) const;
};

#include "Weyl4Cartoon.impl.hpp"

#endif /* WEYL4_HPP_ */
