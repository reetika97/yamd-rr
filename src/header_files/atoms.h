//
// Created by raor on 19.07.24.
//

#ifndef ATOMS_H
#define ATOMS_H

#include "types.h"

class Atoms {
  public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Masses_t masses;
    Energies_t energies;

    Atoms(const int nb_atoms) :
          positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms},
          masses{nb_atoms},energies{nb_atoms}{
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setConstant(1);
        energies.setZero();
    }

    Atoms(const Positions_t &p) :
          positions{p}, velocities{3, p.cols()}, forces{3, p.cols()},
          masses{p.cols()},energies{p.cols()}{
        velocities.setZero();
        forces.setZero();
        masses.setConstant(1);
        energies.setZero();
    }

    Atoms(const Positions_t &p, const Velocities_t &v) :
          positions{p}, velocities{v}, forces{3, p.cols()},
          masses{p.cols()},energies{p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.setConstant(1);
        energies.setZero();
    }

    int nb_atoms() const {
        return positions.cols();
    }

    void resize(int new_size){
        positions.conservativeResize(Eigen::NoChange, new_size);
        velocities.conservativeResize(Eigen::NoChange, new_size);
        forces.conservativeResize(Eigen::NoChange, new_size);
        masses.conservativeResize(new_size);
        energies.conservativeResize(new_size);
    }
};

#endif // ATOMS_H
