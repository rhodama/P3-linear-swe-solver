#!/bin/bash
#SBATCH --account=m4776
#SBATCH -C gpu
#SBATCH --gpus 1
#SBATCH --qos=shared
#SBATCH --time=00:30:00
#SBATCH -N 1
#SBATCH -n 1



make gpu
srun ./build/gpu --output gpu.out --nx 10000 --ny 10000 --num_iter 10000 --save_iter 100
