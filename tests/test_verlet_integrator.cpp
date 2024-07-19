//
// Created by raor on 19.07.24.
//

#include "header_files/verlet.h"
#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>  // For rand() and srand()
#include <ctime>    // For time()
#include <cmath>

// Function to generate a random float
float generateRandomFloat() {
    return static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) ;
}

// Function to calculate analytical values of velocity and position
void calculate_analytical(double &r, double &v, double x,double y,double z,double vx,
                          double vy,double vz,double fx,double fy,double fz,double mass){

    double t = 1e-12;
    double ax, ay, az;
    ax = fx/mass; ay = fy/mass; az = fz/mass;

    x = x + vx*t + 0.5*ax*t*t;
    y = y + vy*t + 0.5*ay*t*t;
    z = z + vz*t + 0.5*az*t*t;

    vx = vx + ax*t;
    vy = vy + ay*t;
    vz = vz + az*t;

    v = std::pow((vx*vx+vy*vy+vz*vz), 0.5);

    r = std::pow((x*x+y*y+z*z), 0.5);

}

TEST(VerletTest, RandomTest) {

    double x,y,z, vx, vy, vz, fx, fy, fz, timestep, mass;
    int nb_steps = 1000;
    x = generateRandomFloat(); y= generateRandomFloat(); z= generateRandomFloat();
    vx = generateRandomFloat(); vy= generateRandomFloat(); vz= generateRandomFloat();
    fx = generateRandomFloat(); fy= generateRandomFloat(); fz= generateRandomFloat();
    timestep = 1e-15; mass = 3.27;

    double r_analytical, v_analytical;
    calculate_analytical(r_analytical, v_analytical, x, y, z, vx, vy, vz, fx, fy, fz, mass);

    single_atom_verlet(x, y, z, vx, vy, vz,
                       fx, fy, fz, timestep, mass, nb_steps);

    double r_verlet = std::pow((x*x+y*y+z*z), 0.5);
    double v_verlet = std::pow((vx*vx+vy*vy+vz*vz), 0.5);

    EXPECT_NEAR(r_verlet, r_analytical, 1e-6);
    EXPECT_NEAR(v_verlet, v_analytical, 1e-6);

}
