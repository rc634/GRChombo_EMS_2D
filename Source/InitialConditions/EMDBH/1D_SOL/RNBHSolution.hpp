
#ifndef RNBHSOLUTION_HPP_
#define RNBHSOLUTION_HPP_

// a lightweight and basic polar areal reissner nordstrom black hole
class RNBHSolution
{

  private:

    double G, L, dx, M, Q;
    int gridsize;

    std::vector<double> Ftr;            // F_tr
    std::vector<double> At;            // electric potential
    std::vector<double> omega;          // lapse
    std::vector<double> psi;            // conformal factor

  public:

    RNBHSolution();
    void set_initialcondition_params(EMSBH_params_t m_params_EMSBH,
                        CouplingFunction::params_t m_params_coupling_function,
                                                          const double max_r);
    double get_Ftr_interp(const double r) const;
    double get_At_interp(const double r) const;
    double get_lapse_interp(const double r) const;
    double get_dlapse_interp(const double r) const;
    double get_psi_interp(const double r) const;
    double get_dpsi_interp(const double r) const;
    double calc_Ftr(const double r) const;

    void main();
    void main_polar_areal();
    void main_iso_sc();
};

#include "RNBHSolution.impl.hpp"

#endif /* RNBHSOLUTION_HPP_ */
