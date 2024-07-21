//
// Created by raor on 20.07.24.
//

#ifndef YAMD_RR_FUNCTIONS_H
#define YAMD_RR_FUNCTIONS_H

void energy_conservation_simulation(double sim_length=100, double timestep=0.001,
                         double sigma=1.0,double mass=1.0, double epsilon=1.0);

double berendsen_thermostat_simulation(int nb_atoms= 100, double target_temp=0.3);
double berendsen_execution_time(int nb_atoms= 100, double target_temp=0.3);

#endif // YAMD_RR_FUNCTIONS_H
