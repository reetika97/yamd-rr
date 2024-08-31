#!/bin/bash -x
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --time=04:00:00
#SBATCH -J YAMD
#SBATCH --mem=6gb
#SBATCH --export=ALL
#SBATCH --partition=multiple

module load compiler/gnu/13.3
module load mpi/openmpi/5.0

echo "Running on ${SLURM_JOB_NUM_NODES} nodes with ${SLURM_JOB_CPUS_PER_NODE} cores each."
echo "Each node has ${SLURM_MEM_PER_NODE} of memory allocated to this job."
time mpirun ./src_main gold_nano_wire "whisker_6950.xyz" 10e-5 1
time mpirun ./src_main gold_nano_wire "whisker_3120.xyz" 10e-5 1
time mpirun ./src_main gold_nano_wire "whisker_840.xyz" 10e-5 1

