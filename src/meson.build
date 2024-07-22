lib_sources = [  # All source files (excluding headers)
    'hello.cpp', 'verlet.cpp', 'lj_direct_summation.cpp', 'berendsen.cpp', 'neighbors.cpp'
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
    'neighbors.cpp', 'equilibration_with_rc.cpp',
    include_directories : [lib_incdirs],
    link_with : [lib],
    dependencies : [eigen, mpi]
)

fs = import('fs')
fs.copyfile('data_files/lj54.xyz')