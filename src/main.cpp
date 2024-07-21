//
// Created by raor on 19.07.24.
//
#include "header_files/xyz.h"
#include "header_files/functions.h"
#include <iostream>

int main(int argc, char* argv[]) {

    std::string pgm_selection = "berendsen_simulation"; //"Default";
    /*if(argc==1) {
        std::cout << "Select_program: " && std::cin >> pgm_selection;
        std::cout << pgm_selection << " is being executed." << std::endl;
    }
    else if(argc==2 ){
        pgm_selection = argv[1];
        std::cout <<pgm_selection<< " is being executed." << std::endl;
    }*/
    if(pgm_selection=="energy_conservation_simulation"){

        energy_conservation_simulation();

    }
    else if(pgm_selection=="berendsen_simulation"){

        auto time_sim=berendsen_thermostat_simulation(50);
        std::cout << "Execution time: " << time_sim << " seconds" << std::endl;

    }

    else {
        std::cout << "Invalid argument(s)" << std::endl;
        exit(0);
    }

    return 0;
}