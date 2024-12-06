/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EMSBHSOLUTION_READ_HPP_)
#error "This file should only be included through EMSBHSolution_read.hpp"
#endif

#ifndef EMSBHSOLUTION_READ_IMPL_HPP_
#define EMSBHSOLUTION_READ_IMPL_HPP_

EMSBHSolution_read::EMSBHSolution_read()
{
}

void EMSBHSolution_read::main(std::string a_data_path)
{
    m_data_path = a_data_path;
    read_from_file();

    gridpoints = r.size();
    L = r[gridpoints-1]-r[0];
    dx = L/( (double) (gridpoints-1));


}

void EMSBHSolution_read::read_from_file()
{
    std::ifstream readfile;

    // maybe guard this ??
    readfile.open(m_data_path);



    std::string line, word;
    while (getline(readfile,line))
    {
        if (line[0]!='#')
        {
            std::vector<int> indices;
            indices.push_back(0);
            for (int i=0; i<line.length(); i++)
            {
                if (line[i]=='\t')
                {
                    indices.push_back(i);
                }
            }

            if (indices.size() < 20)
            {
                pout() << "Warning, cannot properly read initial data dat!"
                       << std::endl
                       << "Number of columns found : " << indices.size()
                       << std::endl;
            }
            else
            {
                r_dummy = extract_string_from_line(line,indices[0],indices[1]);
                phi_dummy = extract_string_from_line(line,indices[1],indices[2]);
                pi_dummy = extract_string_from_line(line,indices[4],indices[5]);
                Er_dummy = extract_string_from_line(line,indices[7],indices[8]);
                a_r_dummy = extract_string_from_line(line,indices[8],indices[9]);
                varphi_dummy = extract_string_from_line(line,indices[9],indices[10]);
                a_dummy = extract_string_from_line(line,indices[10],indices[11]);
                b_dummy = extract_string_from_line(line,indices[11],indices[12]);
                lapse_dummy = extract_string_from_line(line,indices[12],indices[13]);
                shift_dummy = extract_string_from_line(line,indices[13],indices[14]);
                X_dummy = extract_string_from_line(line,indices[14],indices[15]);
                K_dummy = extract_string_from_line(line,indices[15],indices[16]);
                Aa_dummy = extract_string_from_line(line,indices[16],indices[17]);
                Br_dummy = extract_string_from_line(line,indices[18],indices[19]);

                insert_strings_to_vecs();
            }
        }
    }
    readfile.close();
}

void EMSBHSolution_read::insert_strings_to_vecs()
{
    r.push_back(std::stod(r_dummy));
    phi.push_back(std::stod(phi_dummy));
    pi.push_back(std::stod(pi_dummy));
    Er.push_back(std::stod(Er_dummy));
    a_r.push_back(std::stod(a_r_dummy));
    varphi.push_back(std::stod(varphi_dummy));
    a.push_back(std::stod(a_dummy));
    b.push_back(std::stod(b_dummy));
    lapse.push_back(std::stod(lapse_dummy));
    shift.push_back(std::stod(shift_dummy));
    X.push_back(std::stod(X_dummy));
    K.push_back(std::stod(K_dummy));
    Aa.push_back(std::stod(Aa_dummy));
    Br.push_back(std::stod(Br_dummy));
}

std::string EMSBHSolution_read::extract_string_from_line(std::string input,
                                                          int imin, int imax)
{
    std::string dummy = "";
    for (int i=imin; i<imax; i++)
    {
        dummy += input[i];
    }
    return dummy;
}



// 4th order error (cubic interpolation) for field. shouts if asked to fetch a
// value outside the ode solution
double EMSBHSolution_read::get_value_interp(const std::vector<double>& in,
                                                          const double r_) const
{
    // index of 2nd (out of 4) gridpoints used for interpolation
    int iter = (int)floor( (r_-r[0]) / dx);

    // fraction from midpoint of two values,
    // a = +- 1/2 is the nearest gridpoints
    double a = ((r_-r[0]) / dx) - floor((r_-r[0]) / dx) - 0.5;

    double interpolated_value = 0, f1, f2, f3, f4;

    // conditionl/ternary imposing zero gradeint at r=0
    f1 =((iter == 0) ? in[1] : in[iter - 1]);
    f2 = in[iter];
    f3 = in[iter + 1];
    f4 = in[iter + 2];

    if (iter > gridpoints - 3)
    {
        std::cout << "Requested Value outside initial data domain!"
                  << std::endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}

// 4th order error (cubic interpolation) for field. shouts if asked to fetch a
// value outside the ode solution
double EMSBHSolution_read::get_deriv_interp(const std::vector<double>& in,
                                                          const double r_) const
{
    // index of 2nd (out of 4) gridpoints used for interpolation
    int iter = (int)floor((r_-r[0] )/ dx);

    // fraction from midpoint of two values,
    // a = +- 1/2 is the nearest gridpoints
    double a = (r_ / dx) - floor(r_ / dx) - 0.5;

    double interpolated_value = 0, f1, f2, f3, f4;

    // conditionl/ternary imposing zero gradeint at r=0
    f1 =((iter == 0) ? in[1] : in[iter - 1]);
    f2 = in[iter];
    f3 = in[iter + 1];
    f4 = in[iter + 2];

    if (iter > gridpoints - 3)
    {
        std::cout << "Requested Value outside initial data domain!"
                  << std::endl;
    }

    // do the cubic spline (for gradient now), from mathematica script written
    // by Robin (rc634@cam.ac.uk)
    interpolated_value =
        (1. / (24. * dx)) *
        ((f1 - 27. * f2 + 27. * f3 - f4) + 12. * a * (f1 - f2 - f3 + f4) -
         12. * a * a * (f1 - 3. * f2 + 3. * f3 - f4));
    return interpolated_value;
}



#endif /* EMSBHSOLUTION_READ_IMPL_HPP_ */
