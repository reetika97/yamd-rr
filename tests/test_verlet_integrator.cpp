//
// Created by raor on 19.07.24.
//

#include "header_files/verlet.h"
#include "header_files/types.h"
#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>  // For rand() and srand()
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


// Multi atom VerletTest
TEST(VerletTest, MultiAtomTest) {

    Eigen::ArrayXd vel(3);
    Eigen::ArrayXd pos(3);
    Eigen::ArrayXd f(3);

    pos(0) = generateRandomFloat(); pos(1)= generateRandomFloat(); pos(2)= generateRandomFloat();
    vel(0) = generateRandomFloat(); vel(1)= generateRandomFloat(); vel(2)= generateRandomFloat();
    f(0) = generateRandomFloat(); f(1)= generateRandomFloat(); f(2)= generateRandomFloat();

    int nb_steps = 1000;
    double timestep, mass;
    timestep = 1e-15; mass = 3.27;

    int nb_atoms = 4;
    Positions_t positions(3, nb_atoms);
    Velocities_t velocities(3, nb_atoms);
    Forces_t forces(3, nb_atoms);


    auto vx{velocities.row(0)};
    auto vy{velocities.row(1)};
    auto vz{velocities.row(2)};
    auto fx{forces.row(0)};
    auto fy{forces.row(1)};
    auto fz{forces.row(2)};
    auto x{positions.row(0)};
    auto y{positions.row(1)};
    auto z{positions.row(2)};

    vx = vel(0); vy = vel(1); vz = vel(2);
    x = pos(0); y = pos(1); z = pos(2);
    fx = f(0); fy = f(1); fz = f(2);


    single_atom_verlet(pos(0), pos(1), pos(2), vel(0),
                       vel(1), vel(2), f(0), f(1), f(2),
                       timestep, mass, nb_steps);

    multi_atom_verlet(positions, velocities, forces, timestep, mass, nb_steps);

    EXPECT_NEAR(velocities(0), vel(0), 1e-6);
    EXPECT_NEAR(velocities(1), vel(1), 1e-6);
    EXPECT_NEAR(velocities(2), vel(2), 1e-6);
    EXPECT_NEAR(positions(0), pos(0), 1e-6);
    EXPECT_NEAR(positions(1), pos(1), 1e-6);
    EXPECT_NEAR(positions(2), pos(2), 1e-6);

}
