//
// Created by raor on 22.07.24.
//
#include "header_files/xyz.h"
#include <header_files/berendsen.h>
#include <iostream>
#include <chrono>
#include <header_files/lattice.h>
#include <header_files/lj_direct_summation.h>
#include <header_files/verlet.h>
#include <header_files/types.h>
#define save_interval 1000

/**
 * @brief Simulates the equilibration of molecular dynamics system of atoms
 * using Berendsen thermostat and Lennard-Jones potential with
 * cutoff radius and measures the simulation duration. The simulation runs in two phases:
 *
 * 1. **Equilibration Phase:** The Berendsen thermostat is applied to adjust the system temperature to a target value.
 * 2. **Energy conservation Phase:** The thermostat is turned off, allowing the system to relax naturally, with energy conservation monitored.
 *
 * The function logs potential energy, kinetic energy, total energy at regular intervals, and optionally writes trajectories to files.
 *
 * @param nb_atoms The number of atoms in the simulation. Defaults to 50.
 * @param target_temp The target temperature for the system during the temperature control phase. Defaults to 0.3.
 * @param write_to_file If true, the energy data and atom trajectories are written to files ("E.csv" and "traj2.xyz"). Defaults to true.
 *
 * All other parameters can be changed directly in the code
 *
 * @return double The duration of the simulation in seconds.
 */

double equilibration_with_rc(int nb_atoms=50, double target_temp=0.3,
                             bool write_to_file=true){

    auto start = std::chrono::high_resolution_clock::now();

    double sim_length,timestep, sigma=1, mass=1, epsilon=1, kb=1.0, T, rc=5.0;

    auto atoms = lattice_init(nb_atoms);

    std::cout<<"Simulation begins: "<<nb_atoms<<" atoms"<<std::endl;

    timestep = 0.001;
    sim_length = 100;

    auto Epot = lj_direct_summation_rc(atoms, rc);
    auto Ekin = ekin_direct_summation(atoms, mass);
    double Etot = Epot + Ekin;

    std::ofstream traj,E;
    if(write_to_file){
        traj.open("traj2.xyz");
        E.open("E.csv");
        E << "Epot;Ekin;Etot" << std::endl;
    }

    for(int i=0; i<sim_length/timestep; i++) {

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);
        Epot = lj_direct_summation_rc(atoms, rc);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);

        Ekin = ekin_direct_summation(atoms, mass);

        auto T_old = 2 * Ekin / (atoms.nb_atoms() * kb *3);
        berendsen_thermostat(atoms,T_old,timestep,
                             100 * timestep, target_temp);
        Ekin = ekin_direct_summation(atoms, mass);
        T = 2 * Ekin / (atoms.nb_atoms() * kb *3);

        if(i%save_interval == 0 && write_to_file){
            Etot = Epot + Ekin;
            std::cout<<"T-old: "<<T_old<<std::endl;
            std::cout<<"T-new: "<<T<<std::endl;
            std::cout<<"Energy(pot+kin): "<<Epot<<" + "<<Ekin<<" = "<<Etot<<std::endl;
            E<<Epot<<";"<<Ekin<<";"<<Etot<<std::endl;
            write_xyz(traj, atoms); // All trajectories in one file.
        }
    }
    std::cout<<"Berendsen Thermostat Switched off(Equilibrating). "
                 "Energy is conserved."<<std::endl;

    for(int i=0; i<sim_length/timestep; i++) {

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);
        Epot = lj_direct_summation_rc(atoms, rc);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        Ekin = ekin_direct_summation(atoms, mass);

        if(i%save_interval == 0 && write_to_file){
            Etot = Epot + Ekin;
            std::cout<<"Energy(pot+kin): "<<Epot<<" + "<<Ekin<<" = "<<Etot<<std::endl;
            E<<Epot<<";"<<Ekin<<";"<<Etot<<std::endl;
            write_xyz(traj, atoms); // All trajectories in one file.
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return duration.count();
}