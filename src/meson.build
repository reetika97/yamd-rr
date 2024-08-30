lib_sources = [  # All source files (excluding headers)
    'hello.cpp', 'verlet.cpp', 'lj_direct_summation.cpp', 'berendsen.cpp',
    'neighbors.cpp','ducastelle.cpp'
]

lib = static_library(
    'lib',  # Unique name
    lib_sources,
    dependencies : [eigen, mpi]
)

lib_incdirs = include_directories('.')

executable(
    'src_main',
    'main.cpp', 'xyz.cpp',
    'Energy_conservation_LJunits.cpp', 'lj_direct_summation.cpp','lattice.cpp',
    'berendsen_thermostat_simulation.cpp', 'berendsen.cpp',
    'neighbors.cpp', 'equilibration_with_rc.cpp', 'domain.cpp',
    'gold_melting_point.cpp', 'preheat_atom_cluster.cpp', 'ducastelle.cpp', 'energy_conservation_mpi.cpp',
    'gold_nanowire.cpp',

    include_directories : [lib_incdirs],
    link_with : [lib],
    dependencies : [eigen, mpi]
)

fs = import('fs')
fs.copyfile('data_files/lj54.xyz')
fs.copyfile('data_files/cluster_923.xyz')
fs.copyfile('data_files/cluster_3871.xyz')
fs.copyfile('data_files/whisker_small.xyz')
fs.copyfile('data_files/whisker_large.xyz')
fs.copyfile('data_files/whisker_small_T0.xyz')
fs.copyfile('data_files/icosahedrons/55_traj.xyz')
fs.copyfile('data_files/icosahedrons/147_traj.xyz')
fs.copyfile('data_files/icosahedrons/309_traj.xyz')
fs.copyfile('data_files/icosahedrons/561_traj.xyz')
fs.copyfile('data_files/icosahedrons/1415_traj.xyz')
fs.copyfile('data_files/icosahedrons/2869_traj.xyz')
fs.copyfile('data_files/icosahedrons/6525_traj.xyz')
fs.copyfile('data_files/icosahedrons/10179_traj.xyz')
