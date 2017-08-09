#include "jmt.h"

#include <fstream>
#include <iostream>

int main()
{
    double maneuver_time = 5, dt = 0.02;
    double car_s = 124.834, car_d = 6.16483, car_speed = 0., v = 25.;
    auto coeff_s = jmt({car_s, car_speed, 0}, {car_s + 100, 20, 0}, maneuver_time);
    auto coeff_d = jmt({car_d, 0, 0}, {7, 0, 0}, maneuver_time);
    for(int i = 0; i < maneuver_time / dt; i++)
    {
        std::cout << poly(i * dt, coeff_s) << " " << poly(i * dt, coeff_d) << "\n";
    }

    return 0;
}
