//
// Created by raor on 15.06.24.
//
#include "header_files/atoms.h"
#include <iostream>
#include <cmath>


double berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                          double relaxation_time, double target_temp=0.3){


    double rescale_factor;
    rescale_factor = 1 + (((target_temp/temperature)-1)*(timestep/relaxation_time));
    rescale_factor = std::pow(rescale_factor, 0.5);
    //std::cout<<"lambda: "<<rescale_factor<<std::endl;
    atoms.velocities = atoms.velocities * rescale_factor;

    return rescale_factor;


}

void velocity_rescale_energy(Atoms &atoms, double temperature, double del_temp=25){

    double rescale_factor;
    rescale_factor = 1+(del_temp/temperature);
    rescale_factor = std::pow(rescale_factor, 0.5);
    //std::cout<<"lambda: "<<rescale_factor<<std::endl;
    atoms.velocities = atoms.velocities * rescale_factor;


}