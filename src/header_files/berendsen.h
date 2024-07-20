//
// Created by raor on 15.06.24.
//

#ifndef BRENDNSEN_H
#define BRENDNSEN_H

#include "atoms.h"
void berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                          double relaxation_time, double target_temp=0.3);

void velocity_rescale_energy(Atoms &atoms, double temperature, double del_temp=25);

#endif // BRENDNSEN_H
