//
// Created by raor on 19.07.24.
//

#ifndef VERLET_H
#define VERLET_H

#include <Eigen/Dense>

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep, double mass);
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep, double mass);

void single_atom_verlet(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                        double fx, double fy, double fz, double timestep, double mass, int nb_steps);

void verlet_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  const Eigen::Array3Xd &forces, double timestep, double mass);
void verlet_step2(Eigen::Array3Xd &velocities, const Eigen::Array3Xd &forces,
                  double timestep, double mass);

void multi_atom_verlet(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                       const Eigen::Array3Xd &forces, double timestep, double mass, int nb_steps);

#endif // VERLET_H
