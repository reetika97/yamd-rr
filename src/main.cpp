//
// Created by raor on 19.07.24.
//
#include <iostream>
#include <header_files/xyz.h>
#include <header_files/functions.h>


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

    else if(pgm_selection=="equilibration_rc_execution_time"){

        std::ofstream exec_time("run_time_with_rc.csv");
        exec_time<<"num_atoms;exec_time"<<std::endl;
        for(int i=5; i<=150; i+=5){

            auto time_sim=equilibration_with_rc(i, 0.3, false);
            std::cout << "Execution time: " << time_sim << " seconds" << std::endl;
            exec_time<<i<<";"<<time_sim<<std::endl;

        }

    }

    else if(pgm_selection=="equilibration_with_rc"){

        auto time_sim=equilibration_with_rc(50);
        std::cout << "Execution time: " << time_sim << " seconds" << std::endl;

    }

    else if(pgm_selection=="gold_melting_point"){

        gold_melting_point("cluster_923.xyz", true);
        gold_melting_point("923_heated_cluster.xyz");
        gold_melting_point("cluster_3871.xyz", true);
        gold_melting_point("3871_heated_cluster.xyz");

    }

    else {
        std::cout << "Invalid argument(s)" << std::endl;
        exit(0);
    }

    return 0;
}