#ifndef JMT_H
#define JMT_H

#include "Eigen-3.3/Eigen/Dense"

#include <vector>

std::vector<double> jmt(std::vector<double> start, std::vector <double> end, double T)
{
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT
    an array of length 6, each value corresponding to a coefficent in the polynomial
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
    */

    const auto T2 = T * T;
    const auto T3 = T * T2;
    const auto T4 = T * T3;
    const auto T5 = T * T4;

    auto m = Eigen::MatrixXd(3, 3);
    auto rhs = Eigen::VectorXd(3);

    m <<
        T3, T4, T5,
        3 * T2,  4 * T3,  5 * T4,
        6 * T,  12 * T2, 20 * T3;

    rhs <<
        end[0] - (start[0] + start[1] * T + start[2] * T2 / 2.),
        end[1] - (           start[1]     + start[2] * T),
        end[2] - (                          start[2]);

    Eigen::VectorXd alpha(3);
    alpha = m.inverse() * rhs;

    return std::vector<double>{
        start[0],
        start[1],
        start[2] / 2.,
        alpha[0],
        alpha[1],
        alpha[2]};
}

double poly(double x, const std::vector<double> &coeff)
{
    assert(coeff.size() == 6);
    double c0 = coeff[0], c1 = coeff[1], c2 = coeff[2], c3 = coeff[3], c4 = coeff[4], c5 = coeff[5];
    return ((((c5 * x + c4) * x + c3) * x + c2) * x + c1) * x + c0;
}

double dpoly(double x, const std::vector<double> &coeff)
{
    assert(coeff.size() == 6);
    double c0 = coeff[0], c1 = coeff[1], c2 = coeff[2], c3 = coeff[3], c4 = coeff[4], c5 = coeff[5];
    return (((5 * c5 * x + 4 * c4) * x + 3 * c3) * x + 2 * c2) * x + c1;
}

double d2poly(double x, const std::vector<double> &coeff)
{
    assert(coeff.size() == 6);
    double c0 = coeff[0], c1 = coeff[1], c2 = coeff[2], c3 = coeff[3], c4 = coeff[4], c5 = coeff[5];
    return ((20 * c5 * x + 12 * c4) * x + 6 * c3) * x + 2 * c2;
}


#endif
