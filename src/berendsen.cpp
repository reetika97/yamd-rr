//
// Created by raor on 15.06.24.
//
#include "header_files/atoms.h"
#include <iostream>
#include <cmath>

/**
 * @brief Applies the Berendsen thermostat to rescale the velocities of atoms.
 *
 * Simulates the Berendsen thermostat to rescale the velocities of atoms and thereby the temperature.
 * This function calculates a rescale factor based on the difference between the
 * target temperature and the current temperature and the relaxation time.
 *
 * @param atoms `Atoms` object passed by reference containing the current velocities of the atoms.
 * @param temperature The current temperature of the system.
 * @param timestep The time step of the simulation.
 * @param relaxation_time The relaxation time for the Berendsen thermostat. Controls how quickly the temperature approaches the target value.
 * @param target_temp The desired target temperature for the system. Default = 0.3.
 *
 * @return double The rescale factor applied to the atomic velocities.
 *
 */
double berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                          double relaxation_time, double target_temp=0.3){


    double rescale_factor;
    rescale_factor = 1 + (((target_temp/temperature)-1)*(timestep/relaxation_time));
    rescale_factor = std::pow(rescale_factor, 0.5);
    //std::cout<<"lambda: "<<rescale_factor<<std::endl;
    atoms.velocities = atoms.velocities * rescale_factor;

    return rescale_factor;


}

/**
 * @brief Rescales the temperature of system by a specified amount.
 *
 * Rescales the temperature by rescaling velocities of the system by a specified amount.
 * Used to supply external energy into the system.
 *
 *
 * @param atoms `Atoms` object passed by reference containing the current velocities of the atoms.
 * @param temperature The current temperature of the system.
 * @param del_temp The change in temperature to be applied. Default = 25.
 *
 */

void velocity_rescale_energy(Atoms &atoms, double temperature, double del_temp=25){

    double rescale_factor;
    rescale_factor = 1+(del_temp/temperature);
    rescale_factor = std::pow(rescale_factor, 0.5);
    //std::cout<<"lambda: "<<rescale_factor<<std::endl;
    atoms.velocities = atoms.velocities * rescale_factor;

}