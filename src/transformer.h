#ifndef TRANSFORMER_H
#define TRANSFORMER_H

#include <istream>
#include <sstream>
#include <regex>
#include <string>
#include <cmath>

#include <iostream>

struct xy_t
{
    double x, y;
};
xy_t operator+(const xy_t& lhs, const xy_t& rhs) { return xy_t{lhs.x + rhs.x, lhs.y + rhs.y}; };
xy_t operator-(const xy_t& lhs, const xy_t& rhs) { return xy_t{lhs.x - rhs.x, lhs.y - rhs.y}; };
xy_t operator*(double lhs, const xy_t& rhs) { return xy_t{lhs * rhs.x, lhs * rhs.y}; };
xy_t operator/(const xy_t& lhs, double rhs) { return xy_t{lhs.x / rhs, lhs.y / rhs}; };
double operator%(const xy_t& lhs, const xy_t& rhs) {
    auto dx = lhs.x - rhs.x, dy = lhs.y - rhs.y;
    return sqrt(dx * dx + dy * dy);
};

// Cubic Hermite spline with finite difference interpolation https://en.wikipedia.org/wiki/Cubic_Hermite_spline
// p(t) = (2t^3 - 3t^2+1)p0 + (t^3-2t^2+t)m0 + (-2t^3+3t^2)p1 + (t^3-t^2)m1
struct transformer_t
{
    transformer_t(std::istream &input)
    {
        for (std::string line; std::getline(input, line); )
        {
            std::stringstream stream(line);
            double x, y;
            stream >> x >> y;
            p_k.push_back({x, y});
        }

        for (auto it = p_k.begin(); it != p_k.end(); ++it)
        {
            auto prev = std::prev(it == p_k.begin() ? p_k.end() : it);
            auto next = std::next(it) == p_k.end() ? p_k.begin() : std::next(it);

            // https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Finite_difference
            auto d_prev = *it - *prev;
            auto d_next = *next - *it;
            double dt_prev = std::sqrt(d_prev.x * d_prev.x + d_prev.y * d_prev.y);
            double dt_next = std::sqrt(d_next.x * d_next.x + d_next.y * d_next.y);
            m_k.push_back((d_prev / dt_prev + d_next / dt_next) / 2);
        }

        // compute arclength https://en.wikipedia.org/wiki/Curve#Length_of_a_curve
        s_k.push_back(0.);
        for (auto prev = p_k.begin(), it = std::next(p_k.begin()); it != p_k.end(); prev = it, it += 1)
        {
            s_k.push_back(s_k.back() + *it % *prev);
        }
        length = s_k.back() + p_k.front() % p_k.back();
    }

    xy_t p_x(std::size_t k, double s)
    {
        std::size_t next_k = k + 1 == s_k.size() ? 0 : k + 1;
        auto ds = s_k[next_k] - s_k[k];
        auto t = (s - s_k[k]) / ds;
        auto h00 =  2. * t * t * t   - 3. * t * t           + 1.;
        auto h10 =       t * t * t   - 2. * t * t    + t       ;
        auto h01 = -2. * t * t * t   + 3. * t * t               ;
        auto h11 =       t * t * t   -      t * t               ;

        return h00 * p_k[k]
            +  h10 * ds * m_k[k]
            +  h01 * p_k[next_k]
            +  h11 * ds * m_k[next_k];
    }

    xy_t operator()(double s)
    {
        s = fmod(s, length);
        auto it = std::prev(std::upper_bound(s_k.begin(), s_k.end(), s < 0 ? s + length : s));
        std::size_t k = std::distance(s_k.begin(), it);
        return p_x(k, s);
    }

    // xy_t dp_t(std::size_t k, double t)
    // {
    //     auto h00 =  6. * t * t   - 6. * t       ;
    //     auto h10 =  3. * t * t   - 4. * t   + 1.;
    //     auto h01 = -6. * t * t   + 6. * t       ;
    //     auto h11 =  3. * t * t   - 2. * t       ;
    // }

    std::vector<xy_t> p_k;
    std::vector<xy_t> m_k;
    std::vector<double> s_k;
    double length;
};

#endif // TRANSFORMER_H
