//
// Created by raor on 18.06.24.
//
#include "header_files/atoms.h"
#include "header_files/types.h"
#include <cmath>

/**
 * @brief Initializes a simple cubic lattice of atoms based on the number of desired atoms.
 * The lattice is constructed by uniformly placing atoms in a 3D grid, where the lattice constant is set to 1.0 sigma.
 *
 * @param nb_atoms The total number of atoms to be placed in the lattice.
 *
 * @return `Atoms` object containing the positions of the initialized atoms in the lattice.
 *
 * All other parameters can be changed directly in the code
 */

Atoms lattice_init(int nb_atoms){

    double sigma = 1.0; // Lattice constant of 1 sigma
    int N ;
    N = std::ceil(std::pow(nb_atoms, (1/3.))); // N = number of atoms in one dimension

   // nb_atoms positions
    Positions_t pos(3, nb_atoms);

    // Initialize the positions
    for(int l = 0; l < N && nb_atoms > 0; l++){
        for(int k=0; k<= l && nb_atoms > 0; k++)
        {
            for(int j=0; j<=l && nb_atoms > 0; j++ ){

                for(int i=0; i<=l && nb_atoms > 0; i++){

                    //Fill only for newly opened index positions
                    if (i==l || j==l || k==l){

                            pos.col(nb_atoms-1) << i*sigma, j*sigma, k*sigma;
                            nb_atoms--;
                        }

                }
            }
        }
    }
    Atoms atoms(pos);
    return atoms;
}