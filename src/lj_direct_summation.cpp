//
// Created by raor on 04.06.24.
//

#include "header_files/atoms.h"
#include <iostream>
#include <cmath>
#include "header_files/neighbors.h"

/**
 * @brief Computes the Lennard-Jones potential and forces acting on atoms.
 *
 * @param atoms `Atoms` object containing atom positions and forces.
 * @param epsilon Energy parameter Default = 1.0
 * @param sigma Length parameter Default = 1.0
 *
 * @return The total potential energy.
 */

double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0){
    //atoms.energies() is not used or updated in this function!
    //atoms.energies() set to default

    double Epot;
    Epot = 0;
    atoms.forces.setZero();

    for(int i=0; i< atoms.nb_atoms(); i++){

        Eigen::Vector3d fij; //force fij

        for(int j=atoms.nb_atoms()-1; j>i; j--){

            Eigen::Vector3d rij{atoms.positions.col(i) - atoms.positions.col(j)};

            double r = rij.norm(); // magnitude of rij
            if(r==0) std::cout<<"Error";
            Epot += std::pow((sigma/r), 12) - std::pow((sigma/r), 6);

            //Force calculation
            fij = 24 * epsilon * rij/std::pow(r,2) * ((2 * std::pow((sigma/r), 12) )-
                                        std::pow((sigma/r), 6)) ;
            atoms.forces.col(i) += fij.array();
            atoms.forces.col(j) -= fij.array(); //because fji is in opposite direction
        }
    }

    Epot = 4 * epsilon * Epot; // factor of 0.5 is not required as all pairs considered only once
    return Epot;

}

/**
 * @brief Computes the Lennard-Jones potential and forces with a cutoff radius.
 *
 * @param atoms `Atoms` object containing atom positions and forces.
 * @param rc Cutoff radius for the potential calculation.
 * @param epsilon Energy parameter Default = 1.0
 * @param sigma Length parameter Default = 1.0
 *
 * @return The total potential energy.
 */

double lj_direct_summation_rc(Atoms &atoms, double rc, double epsilon = 1.0, double sigma = 1.0){

    double Epot;
    Epot = 0;
    atoms.forces.setZero();
    //modify for loops for j>i
    NeighborList neighbor_list;
    neighbor_list.update(atoms, rc);
    auto Erc = std::pow((sigma/rc), 12) - std::pow((sigma/rc), 6);
    for(auto [i,j]:neighbor_list){

        Eigen::Vector3d fij;

        if(i<j){
            Eigen::Vector3d rij{atoms.positions.col(i) - atoms.positions.col(j)};
            double r = rij.norm();
            if(r==0) std::cout<<"Error";
            Epot += (std::pow((sigma/r), 12) - std::pow((sigma/r), 6) - Erc);
            fij = 24 * epsilon * rij/std::pow(r,2) * ((2 * std::pow((sigma/r), 12) )-
                                                         std::pow((sigma/r), 6)) ;
            atoms.forces.col(i) += fij.array();
            atoms.forces.col(j) -= fij.array();
        }
    }

    Epot = 4 * epsilon * Epot;
    return Epot;

}

/**
 * @brief Computes the kinetic energy of the atoms.
 *
 * @param atoms `Atoms` object containing atom velocities.
 * @param mass Mass of each atom.
 * @param len Number of atoms to consider. Default = all atoms.
 *
 * @return The total kinetic energy.
 */

double ekin_direct_summation(Atoms &atoms, double mass, int len=-10){

    if(len==-10) len = atoms.nb_atoms();

    double Ekin =0.0;
    for(int i=0; i<len; i++){

        Eigen::Vector3d vi{atoms.velocities.col(i)};

        double v = vi.norm();

        Ekin += std::pow(v,2);
    }
    Ekin = 0.5 * mass * Ekin;
    return Ekin;
}