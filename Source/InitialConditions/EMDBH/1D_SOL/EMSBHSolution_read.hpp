
#include <fstream>
#include <sstream>
#include <string>


#ifndef EMSBHSOLUTION_READ_HPP_
#define EMSBHSOLUTION_READ_HPP_

// a lightweight and basic polar areal reissner nordstrom black hole
class EMSBHSolution_read
{

public:

    int gridpoints;
    double dx, L;

    std::vector<double> r;
    std::vector<double> phi;
    std::vector<double> pi;
    std::vector<double> Er;
    std::vector<double> a_r;
    std::vector<double> varphi;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> lapse;
    std::vector<double> shift;
    std::vector<double> X;
    std::vector<double> K;
    std::vector<double> Aa;
    std::vector<double> Br;

    std::string r_dummy;
    std::string phi_dummy;
    std::string pi_dummy;
    std::string Er_dummy;
    std::string a_r_dummy;
    std::string varphi_dummy;
    std::string a_dummy;
    std::string b_dummy;
    std::string lapse_dummy;
    std::string shift_dummy;
    std::string X_dummy;
    std::string K_dummy;
    std::string Aa_dummy;
    std::string Br_dummy;

    std::string m_data_path;

    std::string current_line;

    EMSBHSolution_read();

    double get_value_interp(const std::vector<double>& in, const double r_) const;
    double get_deriv_interp(const std::vector<double>& in, const double r_) const;

    void main(std::string a_data_path);
    void read_from_file();
    std::string extract_string_from_line(std::string input, int imin, int imax);
    void insert_strings_to_vecs();
};

#include "EMSBHSolution_read.impl.hpp"

#endif /* EMSBHSOLUTION_READ_HPP_ */
