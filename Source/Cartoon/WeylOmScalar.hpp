/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEYLOMSCALAR_HPP_
#define WEYLOMSCALAR_HPP_

#include "CCZ4CartoonVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "simd.hpp"
#include <array>

//! Struct for the relevant components of the Riemann tensor with indices down
template <class data_t> struct Riemann_t
{
    Tensor<4, data_t> riemann_dddd; //!< Spatial components of the Riemann
                                    //!< tensor along the evolved directions
    Tensor<3, data_t> riemann_0ddd; //!< One time-like and 3 spatial components
    Tensor<2, data_t> riemann_0d0d; //!< Two time-like and 2 spatial components
    Tensor<2, data_t> riemann_dwdw; //!< Spatial components of the Riemann
                                    //!< tensor with two legs along the sphere
    Tensor<1, data_t> riemann_0wdw; //!< One time-like one spatial and two
                                    //!< sphere components of the Riemann tensor
    data_t riemann_0w0w;            //!< Two time-like and two sphere components
    data_t riemann_wuwu; //!< All components along the transverse sphere
};

//! Struct for the null tetrad
template <class data_t> struct Tetrad_cartoon_t
{
    Tensor<1, data_t> m1; //!< the vector u^i
    Tensor<1, data_t> m2; //!< the vector v^i
    Tensor<1, data_t> m3; //!< the vector w^i
    data_t mw;            //!< vector along the extra dimension;
};

//! Struct for the Newman Penrose scalar matrix
template <class data_t> struct NPScalarMatrix_t
{
    data_t Om22; // 22 component
    data_t Om23; // 23 component
    data_t Om33; // 33 component
    data_t Omww; // component along the extra dimensions
};

//!  Calculates the Weyl scalar matrix for spacetimes without matter content
/*!
 This class calculates the Weyl scalar (real) matrix using definitions
 and tetrad in arXiv:1609.01292 following arXiv:1201.4373.
 */

class WeylOmScalar
{
  public:
    // Use the variable definitions containing the needed quantities
    template <class data_t> using Vars = CCZ4CartoonVars::VarsWithGauge<data_t>;
    template <class data_t>
    using Diff2Vars = CCZ4CartoonVars::Diff2VarsNoGauge<data_t>;

    //! Constructor of class Weyl4Scalar
    /*!
     Takes in the centre for the calculation of the tetrads, and grid spacing
     */
    WeylOmScalar(const std::array<double, CH_SPACEDIM> a_center,
                 const double a_dx)
        : m_center(a_center), m_deriv(a_dx)
    {
    }

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
    const FourthOrderDerivatives m_deriv; //!< for calculating derivs of vars

    //! Calculation of Weyl_4 scalar
    template <class data_t>
    NPScalarMatrix_t<data_t>
    compute_Weyl_Om(const Riemann_t<data_t> &riemann, const Vars<data_t> &vars,
                    const Vars<Tensor<1, data_t>> &d1,
                    const Diff2Vars<Tensor<2, data_t>> &d2,
                    const Coordinates<data_t> &coords) const;

    //! Calculation of the tetrads
    template <class data_t>
    Tetrad_cartoon_t<data_t>
    compute_null_tetrad(const Vars<data_t> &vars,
                        const Coordinates<data_t> &coords) const;

    //! Calculation of the relevant components of the Riemann tensor
    template <class data_t>
    Riemann_t<data_t>
    compute_riemann_tensor(const Vars<data_t> &vars,
                           const Vars<Tensor<1, data_t>> &d1,
                           const Diff2Vars<Tensor<2, data_t>> &d2,
                           const Coordinates<data_t> &coords) const;
};

#include "WeylOmScalar.impl.hpp"

#endif /* WEYLOMSCALAR_HPP_ */
