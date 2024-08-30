//
// Created by raor on 19.07.24.
//
#include <iostream>
#include <header_files/xyz.h>
#include <header_files/functions.h>

/**
 * @brief Main function for selecting and executing different simulation programs.
 *
 * This program allows users to select from energy conservation, Berendsen thermostat,
 * and gold nanowire simulations for single type of particles. Depending
 * on the selection, the corresponding simulation function is executed.
 *
 * @param argc Number of command-line arguments.
 *             If argc == 1, the user is asked to choose selection in the simulation program.
 *             If argc == 2, if second argument is provided it is used as the simulation selection.
 * @param argv Array of command-line arguments.
 *             The second element (argv[1]) is used as the simulation program selection if provided.
 *
 * @return int Returns 0 upon successful execution. Exits with code 0 if invalid input is provided.
 *
 * The available simulations are:
 *  - "energy_conservation_simulation": Executes the energy conservation simulation.
 *  - "berendsen_simulation": Runs the Berendsen thermostat simulation for equilibrating MD system.
 *  - "berendsen_execution_time": Measures execution time of Berendsen thermostat simulation for varying number of atoms.
 *  - "equilibration_with_rc": Runs equilibration with cutoff radius.
 *  - "equilibration_rc_execution_time": Measures execution time of the equilibration with RC simulation for varying number of atoms.
 *  - "gold_melting_point": Executes simulations to determine the melting point of gold clusters.
 *  - "energy_conservation_mpi": Runs the MPI version of the energy conservation simulation.
 *  - "gold_nano_wire": Executes the gold nanowire simulation for calculating stress and strain.
 *
 * If an invalid program selection is made, an error message is displayed and the program exits.
 * Recommended to execute from command line for smoother experience.
 */

int main(int argc, char* argv[]) {

    std::string pgm_selection = "Default"; //"Default";

    if(argc==1) {
        std::cout << "Select_program: " && std::cin >> pgm_selection;
        std::cout << pgm_selection << " is being executed." << std::endl;
    }
    else if(argc>=2 ){
        pgm_selection = argv[1];
        std::cout <<pgm_selection<< " is being executed." << std::endl;
    }

    if(pgm_selection=="energy_conservation_simulation"){

        //Default values of parameters
        double sim_length=100, timestep=0.001, sigma=1.0, mass=1.0, epsilon=1.0;

        if (argc > 2) sim_length = std::atof(argv[2]);
        if (argc > 3) timestep = std::atof(argv[3]);
        if (argc > 4) sigma = std::atof(argv[4]);
        if (argc > 5) mass = std::atof(argv[5]);
        if (argc > 6) epsilon = std::atof(argv[6]);

        energy_conservation_simulation(sim_length, timestep, sigma,
                                       mass, epsilon);

        //energy_conservation_simulation(100,0.0001);
        //energy_conservation_simulation(100,0.005);
        //energy_conservation_simulation(100,0.002);
        //energy_conservation_simulation(100,0.01);
        //energy_conservation_simulation(100,0.02);
        //energy_conservation_simulation(100,0.05);


    }
    else if(pgm_selection=="berendsen_simulation"){

        int nb_atoms= 100; double target_temp=0.3;

        if (argc > 2) nb_atoms = std::atoi(argv[2]);
        if (argc > 3) target_temp = std::atof(argv[3]);

        //auto time_sim=berendsen_thermostat_simulation(100, 0.2);
        auto time_sim = berendsen_thermostat_simulation(nb_atoms,
                                                        target_temp);
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

        int nb_atoms=50; double target_temp=0.3;

        if (argc > 2) nb_atoms = std::atoi(argv[2]);
        if (argc > 3) target_temp = std::atof(argv[3]);

        //auto time_sim=equilibration_with_rc(100, 0.2);
        auto time_sim=equilibration_with_rc(nb_atoms, target_temp);
        std::cout << "Execution time: " << time_sim << " seconds" << std::endl;

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

    else if(pgm_selection=="gold_melting_point"){

        std::string filename; bool preheat_cluster=false;

        if (argc > 2) filename = argv[2];
        if (argc > 3) {
            if (strcmp(argv[3], "true") == 0 || strcmp(argv[3], "1") == 0) preheat_cluster = true;
        }

        gold_melting_point(filename, preheat_cluster);
        //gold_melting_point("cluster_923.xyz", true);
        //gold_melting_point("923_heated_cluster.xyz");
        //gold_melting_point("cluster_3871.xyz", true);
        //gold_melting_point("3871_heated_cluster.xyz");

    }

    else if(pgm_selection=="energy_conservation_mpi"){

        energy_conservation_mpi();

    }

    else if(pgm_selection=="gold_nano_wire"){

        std::string filename="whisker_small.xyz";
        double lx=50, ly=50, lz=140.739;

        if (argc > 2) filename = argv[2];
        if (argc > 3) lx = std::atof(argv[3]);
        if (argc > 4) ly = std::atof(argv[4]);
        if (argc > 5) lz = std::atof(argv[5]);

        preheat_atom_cluster(filename, 10e-5);

        gold_nanowire("3050_heated_cluster.xyz", lx, ly, lz);

    }


    else {
        std::cout << "Invalid argument(s)" << std::endl;
        exit(0);
    }

    return 0;
}