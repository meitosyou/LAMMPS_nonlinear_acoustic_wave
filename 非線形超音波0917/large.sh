#!/bin/sh
#SBATCH --account=HEDISINT
#SBATCH --job-name=lammps
#SBATCH --nodes=256
#SBATCH --ntasks=256
#SBATCH --cpus-per-task=80
#SBATCH --partition=L
#SBATCH --time=23:59:00

srun /project/HEDISINT/lammps/lammps/src/lmp_mpi -in in.lammps
