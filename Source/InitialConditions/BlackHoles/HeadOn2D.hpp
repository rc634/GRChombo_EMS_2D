// Gullstrand-Painlevel coordinates
#ifndef BLACKSTRING_HPP_
#define BLACKSTRING_HPP_

#include "CCZ4CartoonVars.hpp"
//
#include "CCZ4Cartoon.hpp"
//
#include "Cell.hpp"
#include "Coordinates.hpp"
// #include "InitialDataTools.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <array>

//! Class which computes the BlackString
class HeadOn2D
{

  protected:
    BoostedBH::params_t m_bh1_params;
    BoostedBH::params_t m_bh2_params;
    double m_dx;

  public:
    HeadOn2D(BoostedBH::params_t a_bh1_params, BoostedBH::params_t a_bh2_params,
             double a_dx)
        : m_bh1_params(a_bh1_params), m_bh2_params(a_bh2_params), m_dx(a_dx)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const;
};

#include "HeadOn2D.impl.hpp"

#endif /* BLACKSTRING_HPP_ */
