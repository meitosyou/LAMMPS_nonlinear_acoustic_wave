#!/bin/sh
#SBATCH --account=unbilled
#SBATCH --job-name=lammps
#SBATCH --nodes=8
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=10
#SBATCH --partition=S-M
#SBATCH --time=02:59:00

srun /project/HEDISINT/lammps/lammps/src/lmp_mpi -in in.lammps -var rand "$RANDOM"
