# Molecular Dynamics Simulation (Gold Atoms)

## Overview

This project contains various simulation programs designed to execute different molecular dynamics (MD) simulations, such as energy conservation, Berendsen thermostat, gold melting point, and gold nanowire simulations. The program allows users to select and run specific simulations via command-line arguments or through an interactive prompt.

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/reetika97/yamd-rr.git
    ```
2. Compile the program:
    ```bash
    cd <your repository>
    
    # Configure and create build directory
    meson setup builddir --buildtype=release
    
    # Compile
    cd builddir
    meson compile
    
    # Run tests
    meson test

    # Run executables
    cd src
    ./src_main

    #OR for parallel mpi programs (only via command line)
    mpirun -n <#cores> ./src_main <selected_program>
    ```

Ensure that the necessary files, header files and libraries are correctly linked during compilation.

## Usage

The program can be executed via the command line. Depending on the number of command-line arguments provided, the program either prompts for input or uses the arguments directly.

### Command-Line Arguments

- `argc == 1`: The program will prompt the user to select the simulation to run.
- `argc == 2`: The second argument specifies the simulation to execute.
- `argc >= 3`: Additional parameters are used by the selected simulation.

### Available Simulations

| Simulation Program               | Command                    | Description                                                                 | Additional Parameters |
|----------------------------------|----------------------------|-----------------------------------------------------------------------------|-----------------------|
| **Energy Conservation Simulation**   | `energy_conservation_simulation` | Simulates an MD system using LJ-potential and the Velocity Verlet algorithm. Demonstrates the energy conservation in the system.                      | `sim_length` (double), `timestep` (double), `sigma` (double), `epsilon` (double), `mass` (double) |
| **Berendsen Thermostat Simulation**  | `berendsen_simulation`      | Applies the Berendsen thermostat to MD simulation and equilibrates the system.                  | `nb_atoms` (int), `target_temp` (double) |
| **Berendsen Execution Time**         | `berendsen_execution_time`  | Measures the execution time for the Berendsen simulation with varying number of atoms. | None |
| **Equilibration with cutoff radius**            | `equilibration_with_rc`     | Runs MD simulation using LJ-potential with a cutoff radius.                                     | `nb_atoms` (int), `target_temp` (double) |
| **Equilibration with Rc Execution Time**  | `equilibration_rc_execution_time` | Measures the execution time of the equilibration with cutoff radius for varying number of atoms.   | None |
| **Gold Melting Point**               | `gold_melting_point`        | Uses Gupta-Ducastelle potential to simulate melting of gold clusters.                   | `filename` (string), `preheat_cluster` (bool), `temp` (double) |
| **Energy Conservation MPI**          | `energy_conservation_mpi`   | Executes MD simulation using MPI. Demonstrates energy conservation                        | None |
| **Gold Nanowire Simulation**         | `gold_nano_wire`            | Simulates stress and strain on a gold nanowire.                              | `filename` (string), `temp` (double), `del_l` (double) |

## Examples

### Example 1: Energy Conservation Simulation
```bash
# sim_length=100, timestep=0.001, sigma=1.0, epsilon=1.0, mass=1.0 
./src-main energy_conservation_simulation 100 0.001 1.0 1.0 1.0
```
### Example 2: EBerendsen Thermostat Simulation
```bash
# nb_atoms=100, target_temp=0.3
./src-main berendsen_simulation 100 0.3
```

### Example 3: Equilibration with Rc
```bash
# nb_atoms=50, target_temp=0.3
./src-main equilibration_with_rc 50 0.3
```
### Example 4: Gold Melting Point
```bash
# filename="cluster_923.xyz", preheat_cluster=true, temp=600
./src-main gold_melting_point cluster_923.xyz true 600
```

### Example 5: Gold Nanowire Simulation
```bash
# filename="whisker_small.xyz", temp=1e-5, del_l=0.3
./src-main gold_nano_wire whisker_small.xyz 1e-5 0.3
```