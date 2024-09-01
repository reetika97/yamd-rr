//
// Created by raor on 20.07.24.
//

#ifndef YAMD_RR_FUNCTIONS_H
#define YAMD_RR_FUNCTIONS_H

void energy_conservation_simulation(double sim_length=100, double timestep=0.001,
                         double sigma=1.0, double epsilon=1.0, double mass=1.0);

double berendsen_thermostat_simulation(int nb_atoms= 100, double target_temp=0.3,
                                       bool write_to_file=true);

double equilibration_with_rc(int nb_atoms=50, double target_temp=0.3,
                             bool write_to_file=true);

void gold_melting_point(std::string filename);

void energy_conservation_mpi();

void gold_nanowire(std::string filename, double lx, double ly, double lz,
                   double temp, double del_l);

Eigen::Array3d defaults(0.0, 0.0, 0.0);
std::string preheat_atom_cluster(std::string filename, double target_temp = 600,
                                 Eigen::Array3d &domain_lengths = defaults);

#endif // YAMD_RR_FUNCTIONS_H
