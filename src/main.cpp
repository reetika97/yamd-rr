//
// Created by raor on 19.07.24.
//
#include "header_files/xyz.h"
#include "header_files/functions.h"
#include <iostream>

int main(int argc, char* argv[]) {

    std::string pgm_selection = "Default"; //"Default";
    if(argc==1) {
        std::cout << "Select_program: " && std::cin >> pgm_selection;
        std::cout << pgm_selection << " is being executed." << std::endl;
    }
    else if(argc==2 ){
        pgm_selection = argv[1];
        std::cout <<pgm_selection<< " is being executed." << std::endl;
    }
    if(pgm_selection=="energy_conservation_simulation"){

        energy_conservation_simulation();

    }
    else if(pgm_selection=="berendsen_simulation"){

        auto time_sim=berendsen_thermostat_simulation(50);
        std::cout << "Execution time: " << time_sim << " seconds" << std::endl;

    }

    else if(pgm_selection=="berendsen_execution_time"){

        std::ofstream exec_time("run_time_berendsen.csv");
        exec_time<<"num_atoms;exec_time"<<std::endl;
        for(int i=5; i<=150; i+=5){

            auto time_sim=berendsen_thermostat_simulation(i, 0.3, false);
            std::cout << "Execution time: " << time_sim << " seconds" << std::endl;
            exec_time<<i<<";"<<time_sim<<std::endl;

        }

    }

    else if(pgm_selection=="equilibration_with_rc"){

        std::ofstream exec_time("run_time_with_rc.csv");
        exec_time<<"num_atoms;exec_time"<<std::endl;
        for(int i=5; i<=150; i+=5){

            auto time_sim=berendsen_thermostat_simulation(i, 0.3, false);
            std::cout << "Execution time: " << time_sim << " seconds" << std::endl;
            exec_time<<i<<";"<<time_sim<<std::endl;

        }

    }

    else if(pgm_selection=="equilibration_rc_execution_time"){

        auto time_sim=equilibration_with_rc(50, 0.3, false);
        std::cout << "Execution time: " << time_sim << " seconds" << std::endl;

    }

    else {
        std::cout << "Invalid argument(s)" << std::endl;
        exit(0);
    }

    return 0;
}