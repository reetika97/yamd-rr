//
// Created by raor on 04.06.24.
//

#include "header_files/atoms.h"
#include <iostream>
#include <cmath>

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

double ekin_direct_summation(Atoms &atoms, double mass){
    double Ekin =0.0;
    for(int i=0; i<atoms.nb_atoms(); i++){

        Eigen::Vector3d vi{atoms.velocities.col(i)};

        double v = vi.norm(); //velocity magnitude

        Ekin += std::pow(v,2);
    }
    Ekin = 0.5 * mass * Ekin;
    return Ekin;
}