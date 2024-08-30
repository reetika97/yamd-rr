#include "header_files/xyz.h"
#include <cmath>
#include <header_files/atoms.h>
#include <header_files/berendsen.h>
#include <header_files/ducastelle.h>
#include <header_files/lj_direct_summation.h>
#include <header_files/neighbors.h>
#include <header_files/preheat_atom_cluster.h>
#include <header_files/types.h>
#include <header_files/verlet.h>
#include <iostream>
#define RELAX_TIME 1500
#define AVERAGING_WINDOW 500



void gold_melting_point(std::string filename){
    double A = 0.2061, xi = 1.790, p = 10.229, q = 4.036, re = 4.079/sqrt(2);
    double rc = 5.0, timestep = 5, nb_steps = 800000; //1 timestep is 1fs
    double kb = 8.617 * pow(10,-5), mass=196.967*103.6; //eV/K, g/mol
    NeighborList neighbor_list;
    double Epot, Ekin, Etot, T;

    //initialize atom cluster
    auto [names, positions, velocities]
        {read_xyz_with_velocities(filename)};
    Atoms atoms(positions, velocities);
    auto perturbations = atoms.velocities;
    neighbor_list.update(atoms, rc);
    Epot = ducastelle(atoms, neighbor_list, rc, A, xi, p, q, re);
    Ekin = ekin_direct_summation(atoms, mass);
    T = 2 * Ekin / (atoms.nb_atoms() * kb *3);
    std::cout<<atoms.positions<<std::endl;
    std::cout<<atoms.velocities<<std::endl<<std::endl;
    std::cout<<atoms.forces<<std::endl<<std::endl;


    std::cout<<"Simulation begins"<<std::endl;
    std::ofstream traj(std::to_string(atoms.nb_atoms()) +"traj.xyz");
    std::ofstream E(std::to_string(atoms.nb_atoms()) +"E.csv");
    std::ofstream EvsT(std::to_string(atoms.nb_atoms()) +"EvsT.csv");
    E<<"Epot;Ekin;Etot"<<std::endl;
    EvsT<<"Etot;T"<<std::endl;

    int relax_time = RELAX_TIME;
    double avg_E =0.0, avg_T =0.0;
    int x=0;
    for(int i=0; i<nb_steps; i+=timestep) {

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);
        neighbor_list.update(atoms, rc);
        Epot = ducastelle(atoms, neighbor_list, rc, A, xi, p, q, re);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        Ekin = ekin_direct_summation(atoms, mass);
        T = 2 * Ekin / (atoms.nb_atoms() * kb *3);
        Etot = Epot + Ekin;

        //Saves Energy and Trajectory every few steps
        if(i%1000 == 0){
            std::cout<<i<<" T: "<<T<<std::endl;
            std::cout<<"Energy(pot+kin): "<<Epot<<" + "<<Ekin<<" = "<<Etot<<std::endl;
            E<<Epot<<";"<<Ekin<<";"<<Etot<<std::endl;
            write_xyz(traj, atoms); // All trajectories in one file.
        }

        if (relax_time==RELAX_TIME){
            std::cout<<"velocity rescaled"<<std::endl;
            velocity_rescale_energy(atoms, T, 15); //ToDo: Generalize code
            atoms.velocities = atoms.velocities + (perturbations.setRandom()*1e-4);

        }

        if(relax_time<=AVERAGING_WINDOW) {
            x++;
            avg_T += T;
            avg_E += Etot;
        }

        relax_time-=1;
        if(relax_time == 0){

            EvsT<<(avg_E/AVERAGING_WINDOW)<<";"<<(avg_T/AVERAGING_WINDOW)<<std::endl;

            avg_T = 0.0;
            avg_E = 0.0;
            relax_time = RELAX_TIME;

        }
    }

}
