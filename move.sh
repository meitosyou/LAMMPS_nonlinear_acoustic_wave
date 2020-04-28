#!/bin/sh
#SBATCH --account=HEDISINT
#SBATCH --job-name=lammps
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=dev
#SBATCH --time=01:00:00

srun /project/HEDISINT/lammps/lammps/src/lmp_mpi -in in.lammps
