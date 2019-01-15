#!/bin/sh

#SBATCH -t 0:10:00
#SBATCH -J test
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16000

export PATH=/cstor/stanford/toddmtz/users/kgjohn/legion/language:$PATH
module load GCC/4.9.2-binutils-2.25
module load CUDA/8.0.61
module load cuDNN/7.0-CUDA-8.0.61
module load OpenMPI/1.8.5
module load Python/2.7.9
export CONDUIT=ibv
export CC=gcc
export CXX=g++
export GASNET=/cstor/stanford/toddmtz/users/kgjohn/gasnet/release

# run on two nodes
#srun -n 2 -N 2 --tasks-per-node 1 --cpu_bind none regent.py eri.rg -i geom/3x3x3.xyz -b basis/sto-3g.gbs -p 2 -ll:cpu 1 -ll:gpu 1 -ll:csize 14000 -ll:fsize 10000

#run on one node
srun regent.py eri.rg -i geom/4x4x4.xyz -b basis/sto-3g.gbs -o 4x4x4-gpu2.out -p 1 -ll:cpu 1 -ll:gpu 1 -ll:csize 14000 -ll:fsize 10000 

# run with output (NOTE --> output files will be huge for 5x5x5 and above)
#srun regent.py eri.rg -i geom/6x6x6.xyz -b basis/sto-3g.gbs -o 6x6x6.out -p 1 -ll:cpu 1 -ll:gpu 1 -ll:csize 14000 -ll:fsize 10000 
