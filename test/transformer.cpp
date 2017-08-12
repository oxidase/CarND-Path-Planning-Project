#include "transformer.h"
#include "spline.h"
#include "jmt.h"

#include <fstream>
#include <iostream>

int main()
{
    std::ifstream stream("../data/highway_map.csv");
    transformer_t transformer(stream);

    auto x = transformer(6914.15, 3);
    std::cout << x.x << " " << x.y << "\n";

    double maneuver_time = 2;
    auto coeff_s = jmt({6941.8, 20.117, 1.8705e-06}, {6982, 20.117, 0}, maneuver_time);
    auto coeff_d = jmt({10, 0, 0}, {10, 0, 0}, maneuver_time);

    double s = poly(maneuver_time, coeff_s);
    double d = poly(maneuver_time, coeff_d);
    auto pt = transformer(s, d);


    std::vector<double> t_path, x_path, y_path;
    t_path.push_back(0);
    x_path.push_back(779.8063);
    y_path.push_back(1125.686);

    t_path.push_back(0.02);
    x_path.push_back(780.204);
    y_path.push_back(1125.681);

    t_path.push_back(0.08);
    x_path.push_back(781.4102);
    y_path.push_back(1125.664);

    for (double step = maneuver_time / 8, t = maneuver_time / 2; t <= maneuver_time; t += step)
    {
        double s = poly(t, coeff_s);
        double d = poly(t, coeff_d);
        auto pt = transformer(s, d);
        t_path.push_back(t);
        x_path.push_back(pt.x);
        y_path.push_back(pt.y);
        std::cout << t_path.back() << " " << x_path.back() << " " << y_path.back() << " " << s << " " << d << "\n";
    }

    // Interpolate trajectory
    std::vector<double> x_traj, y_traj;
    tk::spline x_spline, y_spline;
    x_spline.set_points(t_path, x_path);
    y_spline.set_points(t_path, y_path);
    for(double t = 0; t <= maneuver_time; t += 0.02)
    {
        x_traj.push_back(x_spline(t));
        y_traj.push_back(y_spline(t));
    }


    std::cout << x_traj.front() << " " << y_traj.front() << "\n";
    for (std::size_t i=1; i < x_traj.size(); ++i)
    {
        auto speed = sqrt(pow(x_traj[i]-x_traj[i-1],2) + pow(y_traj[i]-y_traj[i-1],2)) / 0.02;
        std::cout << i* 0.02 << "      " << x_traj[i] << " " << y_traj[i] << "      " << speed << "        "  << speed *2.236936   << "\n";
    }

    return 0;
}
