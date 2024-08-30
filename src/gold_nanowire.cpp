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
#define SCALE_INTERVAL 200

void gold_nanowire(std::string filename, double lx, double ly, double lz){

    MPI_Init(nullptr, nullptr);
    auto start = std::chrono::high_resolution_clock::now();

    std::cout<<MPI::comm_size(MPI_COMM_WORLD)<<std::endl;

    //double lz = 144.25, lz0 = 144.25;
    const double lx0=lx, ly0=ly, lz0 = lz;
    Eigen::Array3d domain_length;

    Domain domain(MPI_COMM_WORLD, {lx0, ly0, lz0},
                  {1, 1, MPI::comm_size(MPI_COMM_WORLD)},
                  {0, 0, 1});

    //ToCheck: Why am I getting different values with periodicity 001?

    double A = 0.2061, xi = 1.790, p = 10.229, q = 4.036, re = 4.079/sqrt(2);
    double rc = 7.0, timestep = 5, nb_steps = 20000; //1 timestep is 1fs
    double kb = 8.617 * pow(10,-5), mass=196.967*103.6, strain; //eV/K, g/mol

    double Etot, Epot_loc, Epot_total, Ekin_loc, Ekin_total, ghost_force_g;

    auto [names, positions]
        {read_xyz(filename)};
    Atoms atoms(positions);
    atoms.masses.setConstant(mass);

    std::ofstream traj("traj2.xyz");
    std::ofstream svf("svf.csv");
    svf<<"strain;force_lg\n";

    domain.enable(atoms);

    NeighborList neighbor_list;
    neighbor_list.update(atoms, rc);

    double ghost_force=0.0;
    ducastelle_2(atoms, neighbor_list, domain.nb_local(), ghost_force, rc, A, xi, p, q, re);

    for(int i=0; i<nb_steps; i++) {

        if (i % SCALE_INTERVAL == 0 and i != 0){
            std::cout<<"Scaling at "<<i<<std::endl;
            lz += 0.3;
            domain_length<<50,50,lz;
            domain.scale(atoms, domain_length);
        }
        ghost_force=0;
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep,
                     mass);

        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * rc);
        neighbor_list.update(atoms, rc);

        ducastelle_2(atoms, neighbor_list, domain.nb_local(), ghost_force, rc, A, xi, p, q, re);

        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);

        Epot_loc = atoms.energies.head(domain.nb_local()).sum();
        Ekin_loc= ekin_direct_summation(atoms, mass, domain.nb_local());

        auto T = 2 * Ekin_loc / (atoms.nb_atoms() * kb * 3);
        berendsen_thermostat(atoms,T,timestep,
                             100*timestep, 10e-5);

        MPI_Allreduce(&Epot_loc, &Epot_total, 1,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Ekin_loc, &Ekin_total, 1,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&ghost_force, &ghost_force_g, 1,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


        if ((i+1) % SCALE_INTERVAL == 0) { //SCALE_INTERVAL is same as RELAX TIME
            domain.disable(atoms);
            if(domain.rank() == 0) {
                std::cout<<"Saving at "<<i<<std::endl;
                Etot = Epot_total + Ekin_total;
                write_xyz(traj, atoms);
                strain = (lz-lz0) / lz0;
                std::cout << "Energy(pot+kin): " << Epot_total << " + "
                          << Ekin_total << " = " << Etot << std::endl;
                std::cout << "Strain: " << strain << std::endl;
                std::cout << "F_lg: " << ghost_force_g/MPI::comm_size(MPI_COMM_WORLD) << std::endl;
                svf<< strain <<  ";" << ghost_force_g/MPI::comm_size(MPI_COMM_WORLD) << std::endl;

            }

            domain.enable(atoms);
            domain.update_ghosts(atoms, 2 * rc);
            neighbor_list.update(atoms, rc);
            ducastelle_2(atoms, neighbor_list, domain.nb_local(), ghost_force, rc, A, xi, p, q, re);
        }

    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    MPI_Finalize();
}