//
// Created by raor on 19.07.24.
//

#include "header_files/verlet.h"

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