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

/**
 * @brief Preheats an atomic cluster to a target temperature and equilibrates it.
 *
 * This function simulates the heating and equilibration of an atomic cluster
 * using a Berendsen thermostat. The simulation proceeds in two phases:
 * heating to the target temperature and equilibration.
 *
 * @param filename The name of the .XYZ file containing the initial atomic positions.
 * @param target_temp The target temperature to which the atomic cluster should be heated. Defaults to 600 K.
 *
 * @return string Returns filename of preheated cluster
 * All other parameters can be changed directly in the code.
 */

Eigen::Array3d default_array(0.0, 0.0, 0.0);
std::string preheat_atom_cluster(std::string filename, double target_temp = 600,
                                 Eigen::Array3d &domain_lengths = default_array){

    double A = 0.2061, xi = 1.790, p = 10.229, q = 4.036, re = 4.079/sqrt(2);
    double rc = 7.0, timestep = 5, nb_steps = 5000; //1 timestep is 1fs
    double kb = 8.617 * pow(10,-5), mass=196.967*103.6; //eV/K, g/mol
    NeighborList neighbor_list;
    double Epot, Ekin, Etot, T;

    auto [names, positions]
        {read_xyz(filename)};
    Atoms atoms(positions);
    neighbor_list.update(atoms, rc);
    auto traj_file=std::to_string(atoms.nb_atoms())+ "_heated_cluster.xyz";
    std::ofstream traj(traj_file);
    std::ofstream E(std::to_string(atoms.nb_atoms()) + "_E_preheat.csv");
    E<<"Epot;Ekin;Etot"<<std::endl;

    Epot = ducastelle(atoms, neighbor_list, rc, A, xi, p, q, re);

    for(int i=0; i<nb_steps; i++) {

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep,
                     mass);
        neighbor_list.update(atoms, rc);
        Epot = ducastelle(atoms, neighbor_list, rc, A, xi, p, q, re);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        Ekin = ekin_direct_summation(atoms, mass);
        T = 2 * Ekin / (atoms.nb_atoms() * kb * 3);
        berendsen_thermostat(atoms,T,timestep,
                             100*timestep,target_temp);
        Ekin = ekin_direct_summation(atoms, mass);
        T = 2 * Ekin / (atoms.nb_atoms() * kb * 3);
        Etot = Epot + Ekin;

        // Saves Energy and Trajectory every few steps
        if (i % 100 == 0) {
            std::cout << "T: " << T << std::endl;
            std::cout << "Energy(pot+kin): " << Epot << " + " << Ekin << " = "
                      << Etot << std::endl;
            E << Epot << ";" << Ekin << ";" << Etot << std::endl;

        }
    }
    std::cout<<"Equilibrating : "<<std::endl;

    for(int i=0; i<1000; i++) {

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep,
                     mass);
        neighbor_list.update(atoms, rc);
        Epot = ducastelle(atoms, neighbor_list, rc, A, xi, p, q, re);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        Ekin = ekin_direct_summation(atoms, mass);
        T = 2 * Ekin / (atoms.nb_atoms() * kb * 3);
        Etot = Epot + Ekin;

        if (i % 100 == 0) {
            std::cout << "T: " << T << std::endl;
            std::cout << "Energy(pot+kin): " << Epot << " + " << Ekin << " = "
                      << Etot << std::endl;
            E << Epot << ";" << Ekin << ";" << Etot << std::endl;

        }
    }

    std::cout<<"Writing to file: "<<std::endl;
    std::cout << "T: " << T << std::endl;
    std::cout << "Energy(pot+kin): " << Epot << " + " << Ekin << " = "
              << Etot << std::endl;

    domain_lengths = atoms.positions.rowwise().maxCoeff(); // -
                    // atoms.positions.rowwise().minCoeff() ;

    write_xyz(traj, atoms); // Trajectory of heated cluster

    return(traj_file);

}


