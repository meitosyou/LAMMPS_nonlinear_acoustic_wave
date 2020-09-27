#!/bin/sh
#SBATCH --account=unbilled
#SBATCH --job-name=lammps
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=10
#SBATCH --partition=dev
#SBATCH --time=02:59:00

srun /project/HEDISINT/lammps/lammps/src/lmp_mpi -in in.lammps
