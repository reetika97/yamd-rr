#include <gtest/gtest.h>
#include "header_files/berendsen.h"
#include "header_files/lj_direct_summation.h"

TEST(BerendsenTest, immediate_rescale) {
    constexpr int nb_atoms = 10;
    constexpr double target_temp = 0.3;
    constexpr double timestep = 0.001;
    constexpr double relax_time =
        timestep; // to ensure immediate rescaling of temp
    double kb = 1.0, mass = 1.0;

    Atoms atoms(nb_atoms);
    atoms.positions.setRandom(); // random numbers between -1 and 1
    atoms.velocities.setRandom();

    auto Ekin = ekin_direct_summation(atoms, mass);
    auto T_old = 2 * Ekin / (atoms.nb_atoms() * kb * 3);
    berendsen_thermostat(atoms,T_old,timestep,
                         relax_time,target_temp);
    Ekin = ekin_direct_summation(atoms, mass);
    auto T_new = 2 * Ekin / (atoms.nb_atoms() * kb * 3);
    std::cout<<"T_old: "<<T_old<<std::endl;
    std::cout<<"T_new: "<<T_new<<std::endl;
    EXPECT_NEAR(T_new, target_temp, 10e-06);

}

TEST(BerendsenTest, rescale_factor) {
    constexpr int nb_atoms = 10;
    constexpr double target_temp = 0.3;
    constexpr double timestep = 0.001;
    constexpr double relax_time = 100*timestep; //Good Relaxation_time
    double kb = 1.0, mass = 1.0;

    Atoms atoms1(nb_atoms);
    atoms1.positions.setRandom(); // random numbers between -1 and 1
    atoms1.velocities.setConstant(0.8); //Ensuring that T_old is greater than Target_temp

    auto Ekin = ekin_direct_summation(atoms1, mass);
    auto T_old = 2 * Ekin / (atoms1.nb_atoms() * kb * 3);
    auto rescale_factor = berendsen_thermostat(atoms1,T_old,timestep,
                         relax_time,target_temp);
    Ekin = ekin_direct_summation(atoms1, mass);
    auto T_new = 2 * Ekin / (atoms1.nb_atoms() * kb * 3);
    std::cout<<"T_old: "<<T_old<<std::endl;
    std::cout<<"T_new: "<<T_new<<std::endl;
    EXPECT_TRUE(rescale_factor<1);

    Atoms atoms2(nb_atoms);
    atoms2.positions.setRandom();
    atoms2.velocities.setConstant(0.2); //Ensuring that T_old is greater than Target_temp

    Ekin = ekin_direct_summation(atoms2, mass);
    T_old = 2 * Ekin / (atoms2.nb_atoms() * kb * 3);
    rescale_factor = berendsen_thermostat(atoms2,T_old,timestep,
                         relax_time,target_temp);
    Ekin = ekin_direct_summation(atoms1, mass);
    T_new = 2 * Ekin / (atoms2.nb_atoms() * kb * 3);
    std::cout<<"T_old: "<<T_old<<std::endl;
    std::cout<<"T_new: "<<T_new<<std::endl;
    EXPECT_TRUE(rescale_factor>1);

}