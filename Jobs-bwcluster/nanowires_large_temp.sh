#!/bin/bash -x
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=10:00:00
#SBATCH -J YAMD
#SBATCH --mem=6gb
#SBATCH --export=ALL
#SBATCH --partition=multiple

module load compiler/gnu/13.3
module load mpi/openmpi/5.0

echo "Running on ${SLURM_JOB_NUM_NODES} nodes with ${SLURM_JOB_CPUS_PER_NODE} cores each."
echo "Each node has ${SLURM_MEM_PER_NODE} of memory allocated to this job."
time mpirun ./src_main gold_nano_wire "whisker_large.xyz" 50 1
time mpirun ./src_main gold_nano_wire "whisker_large.xyz" 100 1
time mpirun ./src_main gold_nano_wire "whisker_large.xyz" 300 1
time mpirun ./src_main gold_nano_wire "whisker_large.xyz" 400 1

