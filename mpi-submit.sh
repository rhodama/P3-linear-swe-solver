#!/bin/bash
#SBATCH -A m4776               
#SBATCH -C cpu                 
#SBATCH -q debug         
#SBATCH -t 00:30:00           
#SBATCH -N 4                 
#SBATCH --ntasks-per-node=1    
#SBATCH -J mpi_job     
#SBATCH -o %x-%j.out          
#SBATCH -e %x-%j.err          




make mpi
srun -N 4 -n 4 ./build/mpi --nx 1000 --ny 1000 --num_iter 1000 --output mpi.out --save_iter 20
