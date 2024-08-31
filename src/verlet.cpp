//
// Created by raor on 19.07.24.
//

#include "header_files/verlet.h"

/**
 * @brief Performs the predictor step of the velocity Verlet integration for a single atom.
 *
 * @param x Position of the atom in the x-direction.
 * @param y Position of the atom in the y-direction.
 * @param z Position of the atom in the z-direction.
 * @param vx Velocity of the atom in the x-direction.
 * @param vy Velocity of the atom in the y-direction.
 * @param vz Velocity of the atom in the z-direction.
 * @param fx Force acting on the atom in the x-direction.
 * @param fy Force acting on the atom in the y-direction.
 * @param fz Force acting on the atom in the z-direction.
 * @param timestep Time increment for the integration.
 * @param mass Mass of the atom.
 */

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

/**
 * @brief Performs the corrector step of the velocity Verlet integration for a single atom.
 *
 * @param vx Velocity of the atom in the x-direction.
 * @param vy Velocity of the atom in the y-direction.
 * @param vz Velocity of the atom in the z-direction.
 * @param fx Force acting on the atom in the x-direction.
 * @param fy Force acting on the atom in the y-direction.
 * @param fz Force acting on the atom in the z-direction.
 * @param timestep Time increment for the integration.
 * @param mass Mass of the atom.
 */

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep, double mass) {

    // Velocity update at timestep
    vx += 0.5 * fx * timestep / mass;
    vy += 0.5 * fy * timestep / mass;
    vz += 0.5 * fz * timestep / mass;

}

/**
 * @brief Propagates a single atom using the velocity Verlet method.
 *
 * @param x Position of the atom in the x-direction.
 * @param y Position of the atom in the y-direction.
 * @param z Position of the atom in the z-direction.
 * @param vx Velocity of the atom in the x-direction.
 * @param vy Velocity of the atom in the y-direction.
 * @param vz Velocity of the atom in the z-direction.
 * @param fx Force acting on the atom in the x-direction.
 * @param fy Force acting on the atom in the y-direction.
 * @param fz Force acting on the atom in the z-direction.
 * @param timestep Time increment for the integration.
 * @param mass Mass of the atom.
 * @param nb_steps Number of integration steps to perform.
 */

void single_atom_verlet(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                        double fx, double fy, double fz, double timestep, double mass, int nb_steps){
    //Integrator for nb_steps
    for (int i = 0; i < nb_steps; ++i) {

        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep, mass);

        verlet_step2(vx, vy, vz, fx, fy, fz, timestep, mass);
    }
}

/**
 * @brief Performs the predictor step of the velocity Verlet integration for multiple atoms.
 *
 * @param positions Array of atomic positions.
 * @param velocities Array of atomic velocities.
 * @param forces Array of forces acting on the atoms.
 * @param timestep Time increment for the integration.
 * @param mass Mass of each atom.
 */

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

/**
 * @brief Performs the corrector step of the velocity Verlet integration for multiple atoms.
 *
 * @param velocities Array of atomic velocities.
 * @param forces Array of forces acting on the atoms.
 * @param timestep Time increment for the integration.
 * @param mass Mass of each atom.
 */

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

/**
 * @brief Propagates multiple atoms using the velocity Verlet method..
 *
 * @param positions Array of atomic positions.
 * @param velocities Array of atomic velocities.
 * @param forces Array of forces acting on the atoms.
 * @param timestep Time increment for the integration.
 * @param mass Mass of each atom.
 * @param nb_steps Number of integration steps to perform.
 */

void multi_atom_verlet(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                       const Eigen::Array3Xd &forces, double timestep, double mass, int nb_steps){
    //Integrator for nb_steps
    for (int i = 0; i < nb_steps; ++i) {

        verlet_step1(positions, velocities, forces, timestep, mass);

        verlet_step2(velocities, forces, timestep, mass);
    }
}