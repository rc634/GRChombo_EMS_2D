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
class BlackString
{
  public:
    //! Stuct for the params of the Kerr BH
    struct params_t
    {
        double Lhor;
        double pert_amp;
        double pert_freq;
        double regFrac; //!< Regularization fraction of the initial coordinate
                        //!< singularity
        int L;          // this is the length of the string in the x-direction
        int N1;
        int N2;
        std::array<double, CH_SPACEDIM> center;
    };

  protected:
    double m_dx;
    params_t m_params;
    int m_initial_lapse;
    int m_initial_shift;

  public:
    BlackString(params_t a_params, double a_dx, int a_initial_lapse,
                int a_initial_shift)
        : m_dx(a_dx), m_params(a_params), m_initial_lapse(a_initial_lapse),
          m_initial_shift(a_initial_shift)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const;
};

#include "BlackString.impl.hpp"

#endif /* BLACKSTRING_HPP_ */
