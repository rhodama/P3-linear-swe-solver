#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#include "../common/common.hpp"
#include "../common/solver.hpp"
int nx, ny;
double *d_h, *d_u, *d_v;  
double *d_dh, *d_du, *d_dv;  
double *d_dh1, *d_du1, *d_dv1;  
double *d_dh2, *d_du2, *d_dv2;  
int t=0;
double H, g, dx, dy, dt;
#define BLOCK_SIZE 256 
/**
 * This is your initialization function! We pass in h0, u0, and v0, which are
 * your initial height, u velocity, and v velocity fields. You should send these
 * grids to the GPU so you can do work on them there, and also these other fields.
 * Here, length and width are the length and width of the domain, and nx and ny are
 * the number of grid points in the x and y directions. H is the height of the water
 * column, g is the acceleration due to gravity, and dt is the time step size.
 * The rank and num_procs variables are unused here, but you will need them
 * when doing the MPI version.
 */
void init(double *h0, double *u0, double *v0, double length_, double width_, 
          int nx_, int ny_, double H_, double g_, double dt_,
          int rank_, int num_procs_)
{
    // @TODO: your code here
    int size;
    nx = nx_;
    ny = ny_;
    size=nx*ny*sizeof(double);
    cudaMalloc(&d_h, (nx + 1) * (ny + 1) * sizeof(double));
    cudaMalloc(&d_u, (nx + 1) * ny * sizeof(double));
    cudaMalloc(&d_v, nx * (ny + 1) * sizeof(double));
    cudaMalloc(&d_dh,size);
    cudaMalloc(&d_du,size);
    cudaMalloc(&d_dv,size);
    cudaMalloc(&d_dh1,size);
    cudaMalloc(&d_du1,size);
    cudaMalloc(&d_dv1,size);
    cudaMalloc(&d_dh2,size);
    cudaMalloc(&d_du2,size);
    cudaMalloc(&d_dv2,size);

    cudaMemcpy(d_h, h0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_u, u0, (nx + 1) * ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, v0, nx * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);

    H = H_;
    g = g_;

    dx = length_ / nx;
    dy = width_ / ny;

    dt = dt_;

}

__global__ void compute_and_multistep_kernel(double *h, double *u, double *v,
                                    double *dh, double *dh1, double *dh2,
                                    double *du, double *du1, double *du2,
                                    double *dv, double *dv1, double *dv2,
                                    int nx, int ny, double H, double g,
                                    double dx, double dy, double dt,
                                    double a1, double a2, double a3)
{
    __shared__ double s_h[BLOCK_SIZE];
    __shared__ double s_u[BLOCK_SIZE];
    __shared__ double s_v[BLOCK_SIZE];

    int tid = threadIdx.x;
    int global_idx =blockIdx.x*blockDim.x+threadIdx.x;
    int i = global_idx / ny;
    int j = global_idx % ny;

    if(i < nx && j < ny) {
        s_h[tid] = h(i, j);
        s_u[tid] = u(i, j);
        s_v[tid] = v(i, j);

        __syncthreads();

        double local_du_dx = (u(i + 1, j) - s_u[tid]) / dx;
        double local_dv_dy = (v(i, j + 1) - s_v[tid]) / dy;
        double local_dh_dx = (h(i + 1, j) - s_h[tid]) / dx;
        double local_dh_dy = (h(i, j + 1) - s_h[tid]) / dy;

        double current_dh = -H * (local_du_dx + local_dv_dy);
        double current_du = -g * local_dh_dx;
        double current_dv = -g * local_dh_dy;

        h(i, j) += (a1 * current_dh + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
        u(i + 1, j) += (a1 * current_du + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
        v(i, j + 1) += (a1 * current_dv + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;

        __syncthreads();
        
        dh(i, j) = current_dh;
        du(i, j) = current_du;
        dv(i, j) = current_dv;
        __syncthreads();
    }
}

/**
 * This function computes the ghost cells for the horizontal boundaries.
 * This is done by copying the values from the opposite side of the domain.
 */
__global__ void compute_boundaries_kernel(double *h, double *u, double *v,
                                       int nx, int ny)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int j = tid; j < ny; j += stride) {
        h(nx, j) = h(0, j);
        u(0, j) = u(nx, j);
    }

    __syncthreads();

    for (int i = tid; i < nx; i += stride) {
        h(i, ny) = h(i, 0);
        v(i, 0) = v(i, ny);
    }
}



/**
 * This is your step function! Here, you will actually numerically solve the shallow
 * water equations. You should update the h, u, and v fields to be the solution after
 * one time step has passed.
 */
void step()
{
    int total_threads = nx * ny;
    int num_blocks = (total_threads + BLOCK_SIZE - 1) / BLOCK_SIZE;

    int boundary_blocks = (max(nx, ny) + BLOCK_SIZE - 1) / BLOCK_SIZE;
    compute_boundaries_kernel<<<boundary_blocks, BLOCK_SIZE>>>(d_h, d_u, d_v, nx, ny);

    double a1, a2 = 0.0, a3 = 0.0;
    if (t == 0) {
        a1 = 1.0;
    } else if (t == 1) {
        a1 = 3.0 / 2.0;
        a2 = -1.0 / 2.0;
    } else {
        a1 = 23.0 / 12.0;
        a2 = -16.0 / 12.0;
        a3 = 5.0 / 12.0;
    }

    compute_and_multistep_kernel<<<num_blocks, BLOCK_SIZE>>>(d_h, d_u, d_v,
                                                    d_dh, d_dh1, d_dh2,
                                                    d_du, d_du1, d_du2,
                                                    d_dv, d_dv1, d_dv2,
                                                    nx, ny, H, g, dx, dy, dt,
                                                    a1, a2, a3);

    double *tmp;
    tmp = d_dh2; d_dh2 = d_dh1; d_dh1 = d_dh; d_dh = tmp;
    tmp = d_du2; d_du2 = d_du1; d_du1 = d_du; d_du = tmp;
    tmp = d_dv2; d_dv2 = d_dv1; d_dv1 = d_dv; d_dv = tmp;
    
    t++;
}

/**
 * This is your transfer function! You should copy the h field back to the host
 * so that the CPU can check the results of your computation.
 */
void transfer(double *h)
{
    // @TODO: Your code here
    cudaMemcpy(h, d_h, nx * ny * sizeof(double),
                               cudaMemcpyDeviceToHost);
}

/**
 * This is your finalization function! You should free all of the memory that you
 * allocated on the GPU here.
 */
void free_memory()
{
    // @TODO: Your code here
    cudaFree(d_h);
    cudaFree(d_u);
    cudaFree(d_v);
    cudaFree(d_dh);
    cudaFree(d_du);
    cudaFree(d_dv);
    cudaFree(d_dh1);
    cudaFree(d_du1);
    cudaFree(d_dv1);
    cudaFree(d_dh2);
    cudaFree(d_du2);
    cudaFree(d_dv2);
}  