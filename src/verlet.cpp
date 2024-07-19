//
// Created by raor on 19.07.24.
//

#include "header_files/verlet.h"

//Single Atom Velocity Verlet
void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep, double mass) {

    // Velocity update at timestep/2
    vx += 0.5 * fx * timestep / mass;
    vy += 0.5 * fy * timestep / mass;
    vz += 0.5 * fz * timestep / mass;

    // Position update at timestep
    x += vx * timestep;
    y += vy * timestep;
    z += vz * timestep;

}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep, double mass) {

    // Velocity update at timestep
    vx += 0.5 * fx * timestep / mass;
    vy += 0.5 * fy * timestep / mass;
    vz += 0.5 * fz * timestep / mass;

}

void single_atom_verlet(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                        double fx, double fy, double fz, double timestep, double mass, int nb_steps){
    //Integrator for nb_steps
    for (int i = 0; i < nb_steps; ++i) {

        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep, mass);

        verlet_step2(vx, vy, vz, fx, fy, fz, timestep, mass);
    }
}

//Multi Atom Velocity Verlet
void verlet_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  const Eigen::Array3Xd &forces, double timestep, double mass){

    // Velocity update
    auto vx{velocities.row(0)};
    auto vy{velocities.row(1)};
    auto vz{velocities.row(2)};
    auto fx{forces.row(0)};
    auto fy{forces.row(1)};
    auto fz{forces.row(2)};
    auto x{positions.row(0)};
    auto y{positions.row(1)};
    auto z{positions.row(2)};

    vx += 0.5 * fx * timestep / mass;
    vy += 0.5 * fy * timestep / mass;
    vz += 0.5 * fz * timestep / mass;

    // Position update
    x += vx * timestep;
    y += vy * timestep;
    z += vz * timestep;

}

void verlet_step2(Eigen::Array3Xd &velocities, const Eigen::Array3Xd &forces,
                  double timestep, double mass) {

    auto vx{velocities.row(0)};
    auto vy{velocities.row(1)};
    auto vz{velocities.row(2)};
    auto fx{forces.row(0)};
    auto fy{forces.row(1)};
    auto fz{forces.row(2)};

    vx += 0.5 * fx * timestep / mass;
    vy += 0.5 * fy * timestep / mass;
    vz += 0.5 * fz * timestep / mass;
}

void multi_atom_verlet(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                       const Eigen::Array3Xd &forces, double timestep, double mass, int nb_steps){
    //Integrator for nb_steps
    for (int i = 0; i < nb_steps; ++i) {

        verlet_step1(positions, velocities, forces, timestep, mass);

        verlet_step2(velocities, forces, timestep, mass);
    }
}