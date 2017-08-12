#ifndef FSM_H
#define FSM_H

#include "jmt.h"
#include "json.hpp"
#include "transformer.h"
#include "spline.h"

#include <algorithm>
#include <istream>
#include <sstream>
#include <regex>
#include <string>
#include <cmath>

#if defined(USE_GNUPLOT)
#include "gnuplot-iostream.h"
#endif

static inline int to_lane(double d) { return trunc(d / 4.); }
static inline double lane_center(int lane) { return 2. + lane * 4.; }

struct vehicle_t
{
    int id;
    xy_t x, vx;
    xy_t s, vs;
    int lane() const { return to_lane(s.y); }
};
std::ostream& operator<<(std::ostream& s, const vehicle_t& x)
{
    return s << x.id << ": x = (" << x.x << "), v = (" << x.vx << "), s = (" << x.s << "), vs = (" << x.vs << "), lane " << x.lane();
}


struct fsm_t
{
    fsm_t(std::istream &stream)
        : transformer(stream)
        , lane(-1)
    {
    }

    std::pair<std::vector<double>, std::vector<double>> step(nlohmann::json data)
    {
        std::cout << std::setprecision(5);

        // Main car's localization Data
        double car_x = data["x"];
        double car_y = data["y"];
        double car_s = data["s"];
        double car_d = data["d"];
        double car_yaw = data["yaw"];
        double car_speed = data["speed"];

        std::cout << "car: " << car_x << "," << car_y << " " << car_s << " " << car_d << " " << car_speed << " " << car_yaw << "\n";

        // Previous path data given to the Planner
        auto previous_path_x = data["previous_path_x"];
        auto previous_path_y = data["previous_path_y"];

        std::cout << "previous_path_x: " << previous_path_x << "\n";
        std::cout << "previous_path_y: " << previous_path_y << "\n";

        // Previous path's end s and d values
        double end_path_s = data["end_path_s"];
        double end_path_d = data["end_path_d"];

        // Sensor Fusion Data, a list of all other cars on the same side of the road.
        std::vector<vehicle_t> vehicles;
        std::cout << "   vehicles\n";
        for (auto vehicle : data["sensor_fusion"])
        {
            if (vehicle[6] >= 0 && vehicle[6] <= 12)
            {
                xy_t x{vehicle[1], vehicle[2]};
                xy_t vx{vehicle[3], vehicle[4]};
                xy_t s{vehicle[5], vehicle[6]};
                auto vs = transformer(s.x, vx);
                vehicles.push_back({vehicle[0], x, vx, s, vs,});
                std::cout << "   " << vehicles.back() << "\n";
            }
        }

        // Initialize FSM
        if (lane == -1) lane = to_lane(car_d);

        double speed = maximal_speed;
        bool lane_change = false;

        // Check lane and adjust speed or break
        for (auto vehicle : vehicles)
        {
            double distance = vehicle.s.x - car_s;
            distance += distance < 0 ? transformer.length : 0;
            std::cout << "   distance = " << distance << "\n";

            if (car_d - 3 < vehicle.s.y && vehicle.s.y < car_d + 3)
            {
                std::cout << "      check vehicle\n";
                if (distance < safety_distance)
                {
                    std::cout << "        adjusting speed, distance " << distance << " vehicle.vs.x " << vehicle.vs.x << "\n";
                    if (vehicle.vs.x < speed)
                    {
                        speed = vehicle.vs.x;
                        lane_change = true;
                    }
                }
                if (vehicle.lane() == lane && distance < breaking_distance)
                {
                    std::cout << "        breaking, distance " << distance << " vehicle.vs.x " << vehicle.vs.x << "\n";
                    speed = 0;
                    lane_change = false;
                    break;
                }
            }
        }

        if (lane_change)
        {
            if (can_change_to(lane - 1, car_s, speed, vehicles))
                lane -= 1;
            else if (can_change_to(lane + 1, car_s, speed, vehicles))
                lane += 1;
        }

        std::cout << " new speed " << speed << " new lane " << lane << "\n";


        // Interpolate the current (x,v,a) values from the previous path
        double car_vs = 0., car_as = 0., car_vd = 0., car_ad = 0.;
        if (!previous_path_x.empty())
        {
            auto t = maneuver_time - previous_path_x.size() * dt;
            car_s = poly(t, coeff_s);
            car_d = poly(t, coeff_d);
            car_vs = dpoly(t, coeff_s);
            car_as = d2poly(t, coeff_s);
            car_vd = dpoly(t, coeff_d);
            car_ad = d2poly(t, coeff_d);
        }

        // Compute target position and speed
        double acceleration = fmax(fmin((speed - car_vs) / maneuver_time, maximal_acceleration), -maximal_acceleration);
        double target_s = car_s + (car_vs  + acceleration * maneuver_time / 2.) * maneuver_time;
        double target_vs = car_vs + acceleration * maneuver_time;
        double target_d = lane_center(lane);
        double target_vd = (target_d - car_d) / maneuver_time;

        std::cout << "poly s    " << car_s << ", " << car_vs << ", " << car_as << "   ->   " << target_s << ", " << target_vs << ", 0" << "\n";
        std::cout << "poly d    " << car_d << ", " << car_vd << " " << car_ad << "    ->   " << target_d << ", " << target_vd << ", 0" << "\n";

        coeff_s = jmt({car_s, car_vs, car_as}, {target_s, target_vs, 0.}, maneuver_time);
        coeff_d = jmt({car_d, car_vd, car_ad}, {target_d, target_vd, 0.}, maneuver_time);

        return generate_trajectory(previous_path_x, previous_path_y, car_s, car_d);
    }

    bool can_change_to(int lane, double car_s, double speed, const std::vector<vehicle_t>& vehicles)
    {
        if (lane < 0 || lane > 2)
            return false;

        double min_s = car_s, max_s = car_s + speed * maneuver_time, min_d = lane_center(lane) - 1.5, max_d = lane_center(lane) + 1.5;
        std::cout << "can_change_to:   s " << min_s << ".." << max_s << " " << min_d << ".." << max_d << "\n";
        for (auto vehicle : vehicles)
        {
            double vehicle_min_s = vehicle.s.x - 5, vehicle_max_s = vehicle.s.x + vehicle.vs.x * maneuver_time / 2;
            double vehicle_d0 = vehicle.s.y, vehicle_d1 = vehicle.s.y + vehicle.vs.y * maneuver_time;
            double vehicle_min_d = std::min(vehicle_d0 - 1.5, vehicle_d1 - 1.5), vehicle_max_d = std::max(vehicle_d0 + 1.5, vehicle_d1 + 1.5);
            std::cout << "  vehicle s " << vehicle_min_s << ".." << vehicle_max_s << " " << vehicle_min_d << ".." << vehicle_max_d << "\n";
            if (!(max_s < vehicle_min_s || min_s > vehicle_max_s || max_d < vehicle_min_d || min_d > vehicle_max_d))
                return false;
        }
        return true;
    }

    std::pair<std::vector<double>, std::vector<double>> generate_trajectory(const nlohmann::json &previous_path_x,
                                                                            const nlohmann::json &previous_path_y,
                                                                            double car_s, double car_d)
    {
        double mix_time = maneuver_time / 8;
        std::vector<double> t_path, x_path, y_path;
        if (!previous_path_x.empty())
        {
            for (std::size_t i = 0; i < 5 && i < previous_path_x.size(); ++i)
            {
                t_path.push_back(i * dt);
                x_path.push_back(previous_path_x[i]);
                y_path.push_back(previous_path_y[i]);
            }
        }

        std::vector<double> s_plot{0}, d_plot{car_d};
        for (double step = maneuver_time / 8, t = previous_path_x.empty() ? 0 : maneuver_time / 2; t <= maneuver_time; t += step)
        {
            double s = poly(t, coeff_s);
            double d = poly(t, coeff_d);
            s_plot.push_back(s - car_s);
            d_plot.push_back(d);
            auto pt = transformer(s, d);
            t_path.push_back(t);
            x_path.push_back(pt.x);
            y_path.push_back(pt.y);
        }

        std::cout << "trajectory ";
        for (std::size_t i = 0; i < t_path.size(); ++i)
            std::cout << " " << t_path[i] << "," << x_path[i] << "," << y_path[i];
        std::cout << "\n";

#if defined(USE_GNUPLOT)
        gp << "set xrange [0:12]\n";
        gp << "set yrange [-5:50]\n";
        gp << "set arrow from 4, graph 0 to 4, graph 1 nohead\n";
        gp << "set arrow from 8, graph 0 to 8, graph 1 nohead\n";
        gp << "plot '-' using 2:1 with lines title 'path'\n";
        gp.send1d(std::make_tuple(s_plot, d_plot));


        // gp << "plot '-' using 2:3 with lp title 'path'\n";
        // gp.send1d(std::make_tuple(Tpts, Xpts, Ypts));
#endif

        // Interpolate trajectory
        std::vector<double> x_traj, y_traj;
        tk::spline x_spline, y_spline;
        x_spline.set_points(t_path, x_path);
        y_spline.set_points(t_path, y_path);
        for(double t = 0; t <= maneuver_time; t += dt)
        {
            x_traj.push_back(x_spline(t));
            y_traj.push_back(y_spline(t));
        }


        std::cout << "interpolated trajectory\n\n";
        std::cout << "\nx = ["; for (std::size_t i = 0; i < x_traj.size(); ++i) std::cout << " " << x_traj[i] << ",";
        std::cout << "]\ny = ["; for (std::size_t i = 0; i < y_traj.size(); ++i) std::cout << " " << y_traj[i] << ",";
        std::cout << "]\n";

        return std::make_pair(x_traj, y_traj);
    }

    const double safety_distance = 25; // [m]
    const double breaking_distance = 15; // [m]
    const double maneuver_time = 2.; // time to make a maneuver [s]
    const double maximal_speed = 45. / 2.236936; // [m/s]
    const double maximal_acceleration = 8.; // [m/s^2]
    const double dt = 0.02; // time delta [s]
    const transformer_t transformer;

    // State
    int lane;
    std::vector<double> coeff_s, coeff_d;

#if defined(USE_GNUPLOT)
    Gnuplot gp;
#endif
};

#endif
