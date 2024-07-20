//
// Created by raor on 20.07.24.
//
#include <header_files/xyz.h>
#include <header_files/atoms.h>
#include <iostream>
#include <header_files/lj_direct_summation.h>
#include <header_files/verlet.h>
#include <header_files/types.h>
#define save_interval 1000

void Energy_conservation(double sim_length=100, double timestep=0.001,
                         double sigma=1.0,double mass=1.0, double epsilon=1.0)
{

    //Read the file for atom's position and velocity
    auto [names, positions, velocities]
        {read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms(positions, velocities);

    std::cout<<atoms.positions<<std::endl;
    std::cout<<atoms.forces<<std::endl;

    std::cout<<"Simulation begins"<<std::endl;

    //Initial Potential Energy, Kinetic Energy, Total Energy
    double Epot = lj_direct_summation(atoms, epsilon, sigma);
    double Ekin = ekin_direct_summation(atoms, mass);
    double Etot = Epot + Ekin;
    std::cout<<Epot<<" + "<<Ekin<<" = "<<Etot<<std::endl;

    //Initialize file buffer for Total Energy and Atom Trajectory
    std::ofstream traj("traj.xyz");
    std::ofstream Etot_file("Etot_file.csv");

    Etot_file << "Epot;Ekin;Etot"<<std::endl;

    for(int i=0; i<sim_length/timestep; i++) {

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);
        Epot = lj_direct_summation(atoms); //Calculates forces at (t + del t)
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);

        Ekin = ekin_direct_summation(atoms, mass);
        Etot = Epot + Ekin;

        if(i%save_interval == 0){
            std::cout<<Epot<<" + "<<Ekin<<" = "<<Etot<<std::endl;
            Etot_file<<Epot<<";"<<Ekin<<";"<<Etot<<std::endl;
            write_xyz(traj, atoms); // All trajectories in one file.
        }


    }

}
