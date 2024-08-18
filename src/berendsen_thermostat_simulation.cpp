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

double berendsen_thermostat_simulation(int nb_atoms= 100,
                                       double target_temp=0.3,
                                       bool write_to_file=true){

    auto start = std::chrono::high_resolution_clock::now();

    double sim_length,timestep, sigma=1, mass=1, epsilon=1, kb=1, T;

    auto atoms = lattice_init(nb_atoms);


    std::cout<<"Simulation begins: "<<nb_atoms<<" atoms"<<std::endl;

    timestep = 0.001;
    sim_length = 100;

    auto Epot = lj_direct_summation(atoms);
    auto Ekin = ekin_direct_summation(atoms, mass);
    double Etot = Epot + Ekin;


    std::ofstream traj,E;
    if(write_to_file){
        traj.open("traj_5.xyz");
        E.open("E_5.csv");
        E << "Epot;Ekin;Etot" << std::endl;
    }
    for(int i=0; i<sim_length/timestep; i++) {

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);
        Epot = lj_direct_summation(atoms);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);

        Ekin = ekin_direct_summation(atoms, mass);
        auto T_old = 2 * Ekin / (atoms.nb_atoms() * kb *3);

        berendsen_thermostat(atoms,T_old,timestep,
                             100*timestep, target_temp);
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

    std::cout<<"Berendsen Thermostat Switched-off(Equilibrating)."
                 " Energy will be conserved."<<std::endl;

    for(int i=0; i<sim_length/timestep; i++) {

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);
        Epot = lj_direct_summation(atoms);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);

        Ekin = ekin_direct_summation(atoms, mass);
        T = 2 * Ekin / (atoms.nb_atoms() * kb *3);

        if(i%save_interval == 0 && write_to_file){
            Etot = Epot + Ekin;
            std::cout<<"T: "<<T<<std::endl;
            std::cout<<"Energy(pot+kin): "<<Epot<<" + "<<Ekin<<" = "<<Etot<<std::endl;
            E<<Epot<<";"<<Ekin<<";"<<Etot<<std::endl;
            write_xyz(traj, atoms); // All trajectories in one file.
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return duration.count();

}
