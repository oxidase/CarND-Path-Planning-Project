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
    auto coeff_s = jmt({6907.79, 20.1168, 1.9981e-05}, {6948.02, 20.1168, 0}, maneuver_time);
    auto coeff_d = jmt({10, 0, 0}, {10, 0, 0}, maneuver_time);

    double s = poly(maneuver_time, coeff_s);
    double d = poly(maneuver_time, coeff_d);
    auto pt = transformer(s, d);


    std::vector<double> t_path, x_path, y_path;
    t_path.push_back(0);
    x_path.push_back(752.7057);
    y_path.push_back(1134.456);
    t_path.push_back(0.02);
    x_path.push_back(753.009);
    y_path.push_back(1134.426);

    for (double step = maneuver_time / 8, t = step; t <= maneuver_time; t += step)
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

    return 0;
}
