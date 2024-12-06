
#ifndef RNBHSOLUTION_HPP_
#define RNBHSOLUTION_HPP_

// a lightweight and basic polar areal reissner nordstrom black hole
class RNBHSolution
{

  private:

    double G, M, Q;

  public:

    RNBHSolution();
    void set_initialcondition_params(EMSBH_params_t m_params_EMSBH,
                        CouplingFunction::params_t m_params_coupling_function);

    double calc_Er(const double r) const;
    double calc_psi(const double r) const;
};

#include "RNBHSolution.impl.hpp"

#endif /* RNBHSOLUTION_HPP_ */
