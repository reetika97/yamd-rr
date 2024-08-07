//
// Created by raor on 07.08.24.
//
#include "header_files/xyz.h"
#include <header_files/atoms.h>
#include <header_files/berendsen.h>
#include <iostream>
#include <cmath>
#include <header_files/lj_direct_summation.h>
#include <header_files/verlet.h>
#include <header_files/types.h>
#include <header_files/neighbors.h>
#include <header_files/ducastelle.h>
#include <chrono>
#include <header_files/domain.h>

void gold_nanowire(){

    MPI_Init(nullptr, nullptr);
    auto start = std::chrono::high_resolution_clock::now();

    std::cout<<MPI::comm_size(MPI_COMM_WORLD)<<std::endl;


    Domain domain(MPI_COMM_WORLD, {50, 50, 135},
                  {1, 1, MPI::comm_size(MPI_COMM_WORLD)},
                  {0, 0, 1});

    //ToCheck: Why am I getting different values with periodicity 001?

    double A = 0.2061, xi = 1.790, p = 10.229, q = 4.036, re = 4.079/sqrt(2);
    double rc = 7.0, timestep = 5, nb_steps = 10000; //1 timestep is 1fs
    double kb = 8.617 * pow(10,-5), mass=196.967*103.6; //eV/K, g/mol

    double Etot, Epot_loc, Epot_total, Ekin_loc, Ekin_total;

    auto [names, positions]
        {read_xyz("whisker_small.xyz")};
    Atoms atoms(positions);
    atoms.masses.setConstant(mass);

    std::ofstream traj("traj2.xyz");

    domain.enable(atoms);

    NeighborList neighbor_list;
    neighbor_list.update(atoms, rc);

    ducastelle(atoms, neighbor_list, rc, A, xi, p, q, re);

    for(int i=0; i<nb_steps; i+=timestep) {

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep,
                     mass);

        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * rc);
        neighbor_list.update(atoms, rc);

        ducastelle(atoms, neighbor_list, rc, A, xi, p, q, re);

        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);

        Epot_loc = atoms.energies.head(domain.nb_local()).sum();
        Ekin_loc= ekin_direct_summation(atoms, mass, domain.nb_local());

        MPI_Allreduce(&Epot_loc, &Epot_total, 1,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Ekin_loc, &Ekin_total, 1,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


        if (i % 1000 == 0) {
            domain.disable(atoms);
            if(domain.rank() == 0) {
                Etot = Epot_total + Ekin_total;
                write_xyz(traj, atoms);
                std::cout << "Energy(pot+kin): " << Epot_total << " + "
                          << Ekin_total << " = " << Etot << std::endl;
            }

            domain.enable(atoms);
            domain.update_ghosts(atoms, 2 * rc);
            neighbor_list.update(atoms, rc);
            ducastelle(atoms, neighbor_list, rc, A, xi, p, q, re);
        }

    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    MPI_Finalize();
}