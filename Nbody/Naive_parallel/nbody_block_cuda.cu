#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"

#define BLOCK_SIZE 256
#define SOFTENING 1e-9f

typedef struct { float4 *pos, *vel; } BodySystem;

void randomizeBodies(float4 *pos, float4 *vel, int n) {
  srand(42);
  for (int i = 0; i < n; i++) {
    pos[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    pos[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    pos[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    pos[i].w = 0.0f; // Ignored or set to 0
    
    // Initialize velocity (vx, vy, vz), set w to 0
    vel[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    vel[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    vel[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    vel[i].w = 0.0f; // Ignored or set to 0
  }
}

void savePositionsToFile(BodySystem &p, int nBodies, const char* filepath) {
  FILE *file = fopen(filepath, "w");
  if (file == NULL) {
    printf("Error opening file for writing!\n");
    return;
  }

  for (int i = 0; i < nBodies; i++) {
    // Position from float4
    fprintf(file, "Body %d: Position(%.6f, %.6f, %.6f) ", i, p.pos[i].x, p.pos[i].y, p.pos[i].z);
    // Velocity from float4
    fprintf(file, "Velocity(%.6f, %.6f, %.6f)\n", p.vel[i].x, p.vel[i].y, p.vel[i].z);
  }

  fclose(file);
}

__global__
void bodyForce(float4 *p, float4 *v, float dt, int n) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i < n) {
    float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

    for (int tile = 0; tile < gridDim.x; tile++) {
      __shared__ float3 spos[BLOCK_SIZE];
      float4 tpos = p[tile * blockDim.x + threadIdx.x];
      spos[threadIdx.x] = make_float3(tpos.x, tpos.y, tpos.z);
      __syncthreads();

      for (int j = 0; j < BLOCK_SIZE; j++) {
        float dx = spos[j].x - p[i].x;
        float dy = spos[j].y - p[i].y;
        float dz = spos[j].z - p[i].z;
        float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
        float invDist = rsqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

        Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
      }
      __syncthreads();
    }

    v[i].x += dt*Fx; v[i].y += dt*Fy; v[i].z += dt*Fz;
  }
}

int main(const int argc, const char** argv) {
  
  int nBodies = 30000;
  int nIters = 10;  // simulation iterations
  if (argc > 1) 
  {
  if(argc < 3)  
  {
  	printf("usage: <executable> <number of N-bodies> <iterations>\n");
  return 0;
  }
  else if (argc > 3)
  {
  	printf("usage: <executable> <number of N-bodies> <iterations>\n");
  return 0;
  }
  
  else
  {
  
  nBodies = atoi(argv[1]);
  nIters = atoi(argv[2]);
  }
}
  else 
  {
  printf("usage: <executable> <number of N-bodies> <iterations>\n");
  return 0;
	}
  
  const float dt = 0.01f; // time step
  
  
  int bytes = 2*nBodies*sizeof(float4);
  float *buf = (float*)malloc(bytes);
  BodySystem p = { (float4*)buf, ((float4*)buf) + nBodies };

  randomizeBodies(p.pos, p.vel, nBodies); // Init pos / vel data

  savePositionsToFile(p, nBodies, "gpu_start.txt");

  float *d_buf;
  cudaMalloc(&d_buf, bytes);
  BodySystem d_p = { (float4*)d_buf, ((float4*)d_buf) + nBodies };

  int nBlocks = (nBodies + BLOCK_SIZE - 1) / BLOCK_SIZE;
  double totalTime = 0.0; 

  for (int iter = 1; iter <= nIters; iter++) {
    StartTimer();

    cudaMemcpy(d_buf, buf, bytes, cudaMemcpyHostToDevice);
    bodyForce<<<nBlocks, BLOCK_SIZE>>>(d_p.pos, d_p.vel, dt, nBodies);
    cudaMemcpy(buf, d_buf, bytes, cudaMemcpyDeviceToHost);

    for (int i = 0 ; i < nBodies; i++) { // integrate position
      p.pos[i].x += p.vel[i].x*dt;
      p.pos[i].y += p.vel[i].y*dt;
      p.pos[i].z += p.vel[i].z*dt;
    }

    const double tElapsed = GetTimer() / 1000.0;
    if (iter > 1) { // First iter is warm up
      totalTime += tElapsed; 
    }

    printf("Iteration %d: %.5f seconds\n", iter, tElapsed);

  }
  double avgTime = totalTime / (double)(nIters-1); 


  
  printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
  
  savePositionsToFile(p, nBodies, "gpu_end.txt");
  
  free(buf);
  cudaFree(d_buf);
  return 0;
}
