#!/bin/bash
#SBATCH -A m4776              
#SBATCH -C cpu                
#SBATCH -q debug        
#SBATCH -t 00:30:00         
#SBATCH -N 1                 
#SBATCH --ntasks-per-node=1    
#SBATCH -J serial_job     
#SBATCH -o %x-%j.out          
#SBATCH -e %x-%j.err          





srun ./build/serial --nx 1000 --ny 1000 --num_iter 10000 --output serial_1.out --save_iter 100
