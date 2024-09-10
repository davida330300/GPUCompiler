GPU version:

> module purge
> module load nvhpc
> module load cuda/10.2
> nvcc -c -o nbody_block_cuda.o nbody_block_cuda.cu -O3 -I/usr/local/cuda/include -arch=sm_75
> nvcc nbody_block_cuda.o -o nbody_block_cuda -lcudart -L/usr/local/cuda/lib64
> ./nbody_block_cuda 10000 10

CPU OpenMP version:
> module unload intel/compiler/64/16.0.2/2016.2.181
> gcc -c -Wall -std=c99 -c -o nbody.o nbody.c
> gcc nbody.o -o nbody -lm
> ./nbody 10000 10

// optimization
> gcc -c -Wall -std=c99 –O3 -c -o nbody.o nbody.c
> gcc nbody.o -o nbody -lm

Assignment 1:
