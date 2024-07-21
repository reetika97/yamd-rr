//
// Created by raor on 20.07.24.
//
#include "header_files/xyz.h"
#include <header_files/berendsen.h>
#include <header_files/lattice.h>
#include <header_files/lj_direct_summation.h>
#include <header_files/types.h>
#include <header_files/verlet.h>
#include <iostream>
#include <chrono>
#define save_interval 1000

double berendsen_execution_time(int nb_atoms= 100, double target_temp=0.3){

    auto start = std::chrono::high_resolution_clock::now();

    double sim_length,timestep, sigma=1, mass=1, epsilon=1, kb=1, T;

    auto atoms = lattice_init(nb_atoms);

    std::cout<<"Simulating: "<<nb_atoms<<" atoms"<<std::endl;

    timestep = 0.001;
    sim_length = 100;

    auto Epot = lj_direct_summation(atoms);
    auto Ekin = ekin_direct_summation(atoms, mass);
    double Etot = Epot + Ekin;

    for(int i=0; i<sim_length/timestep; i++) {

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);
        Epot = lj_direct_summation(atoms);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);

        Ekin = ekin_direct_summation(atoms, mass);
        auto T_old = 2 * Ekin / (atoms.nb_atoms() * kb *3);

        berendsen_thermostat(atoms,T_old,timestep,
                             500*timestep, target_temp);
        Ekin = ekin_direct_summation(atoms, mass);
        T = 2 * Ekin / (atoms.nb_atoms() * kb *3);

    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    return duration.count();

}
