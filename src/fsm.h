#ifndef FSM_H
#define FSM_H

#include "transformer.h"
#include "json.hpp"

#include <istream>
#include <sstream>
#include <regex>
#include <string>
#include <cmath>

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
    {
    }

    std::pair<std::vector<double>, std::vector<double>> step(nlohmann::json data)
    {
        // Main car's localization Data
        double car_x = data["x"];
        double car_y = data["y"];
        double car_s = data["s"];
        double car_d = data["d"];
        double car_yaw = data["yaw"];
        double car_speed = data["speed"];

        // Previous path data given to the Planner
        auto previous_path_x = data["previous_path_x"];
        auto previous_path_y = data["previous_path_y"];

        // Previous path's end s and d values
        double end_path_s = data["end_path_s"];
        double end_path_d = data["end_path_d"];

        // Sensor Fusion Data, a list of all other cars on the same side of the road.
        std::vector<vehicle_t> vehicles;
        for (auto vehicle : data["sensor_fusion"])
        {
            if (vehicle[6] >= 0 && vehicle[6] <= 12)
            {
                xy_t x{vehicle[1], vehicle[2]};
                xy_t vx{vehicle[3], vehicle[4]};
                xy_t s{vehicle[5], vehicle[6]};
                auto vs = transformer(s.x, vx);
                vehicles.push_back({vehicle[0], x, vx, s, vs,});
            }
        }

        std::cout << "main car: " << ": x = (" << car_x << ", " << car_y << "), v = " << car_speed << ", s = (" << car_s << "," << car_d<< ")\n";
        for (auto vehicle : vehicles)
            std::cout << vehicle << "\n";

        double v = 25; // velocity [m/s]
        // check lane 2 and adjust speed
        for (auto vehicle : vehicles)
        {
            if (vehicle.lane() == 2 && car_s < vehicle.s.x && vehicle.s.x < car_s + 20)
            {
                v = vehicle.vs.x;
            }
        }

        std::vector<xy_t> steps; // steps in sd-space (s=>x, d=>y)

        for(int i = 0; i < 50; i++)
        {
            steps.push_back({car_s + i * v * dt, 10});
        }

        return transform_to_world(steps);
    }

    std::pair<std::vector<double>, std::vector<double>> transform_to_world(std::vector<xy_t> steps)
    {
        std::vector<double> x, y;
        for (auto step : steps)
        {
            auto pt = transformer(step.x, step.y);
            x.push_back(pt.x);
            y.push_back(pt.y);
        }
        return std::make_pair(x, y);
    }

    const double dt = 0.02; // time delta [s]
    transformer_t transformer;
};

#endif
