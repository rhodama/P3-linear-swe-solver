#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <iostream>
#include <immintrin.h>

#include "../common/common.hpp"
#include "../common/solver.hpp"

#define ALIGNMENT 64
// Here we hold the number of cells we have in the x and y directions
int nx, ny;
// Alignment mask array, used to handle remaining elements

alignas(32) static const long long _pd_mask[4] = {-1LL, -1LL, -1LL, -1LL};


// This is where all of our points are. We need to keep track of our active
// height and velocity grids, but also the corresponding derivatives. The reason
// we have 2 copies for each derivative is that our multistep method uses the
// derivative from the last 2 time steps.
double *h, *u, *v, *dh, *du, *dv, *dh1, *du1, *dv1, *dh2, *du2, *dv2;
double H, g, dx, dy, dt;

/**
 * This is your initialization function!
 */
void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    // TODO: Your code here
    // We set the pointers to the arrays that were passed in
    nx = nx_;
    ny = ny_;

    dh = new (std::align_val_t(64)) double[nx * ny]();  
    du = new (std::align_val_t(64)) double[nx * ny]();
    dv = new (std::align_val_t(64)) double[nx * ny]();

    dh1 = new (std::align_val_t(64)) double[nx * ny]();
    du1 = new (std::align_val_t(64)) double[nx * ny]();
    dv1 = new (std::align_val_t(64)) double[nx * ny]();

    dh2 = new (std::align_val_t(64)) double[nx * ny]();
    du2 = new (std::align_val_t(64)) double[nx * ny]();
    dv2 = new (std::align_val_t(64)) double[nx * ny]();

    h = h0;
    u = u0;
    v = v0;


    H = H_;
    g = g_;

    dx = length_ / nx;
    dy = width_ / ny;

    dt = dt_;
}
/**
 * This function computes the derivative of the height field
 * with respect to time. This is done by taking the divergence
 * of the velocity field and multiplying by -H.
 */
void compute_dh()
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
        }
    }
}

/**
 * This function computes the derivative of the x-component of the
 * velocity field with respect to time. This is done by taking the
 * derivative of the height field with respect to x and multiplying
 * by -g.
 */
void compute_du()
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            du(i, j) = -g * dh_dx(i, j);
        }
    }
}


/**
 * This function computes the derivative of the y-component of the
 * velocity field with respect to time. This is done by taking the
 * derivative of the height field with respect to y and multiplying
 * by -g.
 */
void compute_dv()
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            dv(i, j) = -g * dh_dy(i, j);
        }
    }
}

/**
 * This function computes the next time step using a multistep method.
 * The coefficients a1, a2, and a3 are used to determine the weights
 * of the current and previous time steps.
 */
void multistep(double a1, double a2, double a3)
{
    __m256d va1 = _mm256_set1_pd(a1);
    __m256d va2 = _mm256_set1_pd(a2);
    __m256d va3 = _mm256_set1_pd(a3);
    __m256d vdt = _mm256_set1_pd(dt);
for (int i = 0; i < nx; i++)
{
    int j;
    // Main loop, processing 8 elements at a time
    for (j = 0; j <= ny-8; j += 8)
    {
        // Prefetch data for the next iteration
        _mm_prefetch((const char*)&dh(i, j+8), _MM_HINT_T0);
        _mm_prefetch((const char*)&dh1(i, j+8), _MM_HINT_T0);
        _mm_prefetch((const char*)&dh2(i, j+8), _MM_HINT_T0);
        _mm_prefetch((const char*)&du(i, j+8), _MM_HINT_T0);
        _mm_prefetch((const char*)&du1(i, j+8), _MM_HINT_T0);
        _mm_prefetch((const char*)&du2(i, j+8), _MM_HINT_T0);
        _mm_prefetch((const char*)&dv(i, j+8), _MM_HINT_T0);
        _mm_prefetch((const char*)&dv1(i, j+8), _MM_HINT_T0);
        _mm_prefetch((const char*)&dv2(i, j+8), _MM_HINT_T0);

        // Process the h array
        __m256d vdh1 = _mm256_load_pd(&dh(i, j));
        __m256d vdh2 = _mm256_load_pd(&dh(i, j+4));
        __m256d vdh11 = _mm256_load_pd(&dh1(i, j));
        __m256d vdh12 = _mm256_load_pd(&dh1(i, j+4));
        __m256d vdh21 = _mm256_load_pd(&dh2(i, j));
        __m256d vdh22 = _mm256_load_pd(&dh2(i, j+4));
        __m256d vh1 = _mm256_load_pd(&h(i, j));
        __m256d vh2 = _mm256_load_pd(&h(i, j+4));

        __m256d vtemp1 = _mm256_mul_pd(va1, vdh1);
        __m256d vtemp2 = _mm256_mul_pd(va1, vdh2);
        vtemp1 = _mm256_fmadd_pd(va2, vdh11, vtemp1);
        vtemp2 = _mm256_fmadd_pd(va2, vdh12, vtemp2);
        vtemp1 = _mm256_fmadd_pd(va3, vdh21, vtemp1);
        vtemp2 = _mm256_fmadd_pd(va3, vdh22, vtemp2);
        vtemp1 = _mm256_fmadd_pd(vtemp1, vdt, vh1);
        vtemp2 = _mm256_fmadd_pd(vtemp2, vdt, vh2);

        _mm256_storeu_pd(&h(i, j), vtemp1);
        _mm256_storeu_pd(&h(i, j+4), vtemp2);

        // Process the u array
        __m256d vdu1 = _mm256_load_pd(&du(i, j));
        __m256d vdu2 = _mm256_load_pd(&du(i, j+4));
        __m256d vdu11 = _mm256_load_pd(&du1(i, j));
        __m256d vdu12 = _mm256_load_pd(&du1(i, j+4));
        __m256d vdu21 = _mm256_load_pd(&du2(i, j));
        __m256d vdu22 = _mm256_load_pd(&du2(i, j+4));
        __m256d vu1 = _mm256_load_pd(&u(i+1, j));
        __m256d vu2 = _mm256_load_pd(&u(i+1, j+4));

        vtemp1 = _mm256_mul_pd(va1, vdu1);
        vtemp2 = _mm256_mul_pd(va1, vdu2);
        vtemp1 = _mm256_fmadd_pd(va2, vdu11, vtemp1);
        vtemp2 = _mm256_fmadd_pd(va2, vdu12, vtemp2);
        vtemp1 = _mm256_fmadd_pd(va3, vdu21, vtemp1);
        vtemp2 = _mm256_fmadd_pd(va3, vdu22, vtemp2);
        vtemp1 = _mm256_fmadd_pd(vtemp1, vdt, vu1);
        vtemp2 = _mm256_fmadd_pd(vtemp2, vdt, vu2);

        _mm256_storeu_pd(&u(i+1, j), vtemp1);
        _mm256_storeu_pd(&u(i+1, j+4), vtemp2);

        // Process the v array
        __m256d vdv1 = _mm256_load_pd(&dv(i, j));
        __m256d vdv2 = _mm256_load_pd(&dv(i, j+4));
        __m256d vdv11 = _mm256_load_pd(&dv1(i, j));
        __m256d vdv12 = _mm256_load_pd(&dv1(i, j+4));
        __m256d vdv21 = _mm256_load_pd(&dv2(i, j));
        __m256d vdv22 = _mm256_load_pd(&dv2(i, j+4));
        __m256d vv1 = _mm256_load_pd(&v(i, j+1));
        __m256d vv2 = _mm256_load_pd(&v(i, j+5));

        vtemp1 = _mm256_mul_pd(va1, vdv1);
        vtemp2 = _mm256_mul_pd(va1, vdv2);
        vtemp1 = _mm256_fmadd_pd(va2, vdv11, vtemp1);
        vtemp2 = _mm256_fmadd_pd(va2, vdv12, vtemp2);
        vtemp1 = _mm256_fmadd_pd(va3, vdv21, vtemp1);
        vtemp2 = _mm256_fmadd_pd(va3, vdv22, vtemp2);
        vtemp1 = _mm256_fmadd_pd(vtemp1, vdt, vv1);
        vtemp2 = _mm256_fmadd_pd(vtemp2, vdt, vv2);

        _mm256_storeu_pd(&v(i, j+1), vtemp1);
        _mm256_storeu_pd(&v(i, j+5), vtemp2);
    }

    // Process the remaining elements (if any)
    for (; j < ny; j++)
    {
        // Handle single elements using scalar operations
        h(i, j) += dt * (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j));
        u(i+1, j) += dt * (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j));
        v(i, j+1) += dt * (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j));
    }
}

}


/**
 * This function computes the ghost cells for the horizontal boundaries.
 * This is done by copying the values from the opposite side of the domain.
 */
void compute_ghost_horizontal()
{
    for (int j = 0; j < ny; j++)
    {
        h(nx, j) = h(0, j);
    }
}

/**
 * This function computes the ghost cells for the vertical boundaries.
 * This is done by copying the values from the opposite side of the domain.
 */
void compute_ghost_vertical()
{
    for (int i = 0; i < nx; i++)
    {
        h(i, ny) = h(i, 0);
    }
}

/**
 * This function computes the boundaries for the horizontal boundaries.
 * We do this by copying the values from the opposite side of the domain.
 */
void compute_boundaries_horizontal()
{
    for (int j = 0; j < ny; j++)
    {
        u(0, j) = u(nx, j);
    }
}

/**
 * This function computes the boundaries for the vertical boundaries.
 * We do this by copying the values from the opposite side of the domain.
 */
void compute_boundaries_vertical()
{
    for (int i = 0; i < nx; i++)
    {
        v(i, 0) = v(i, ny);
    }
}

/**
 * This function swaps the buffers for the derivatives of our different fields.
 * This is done so that we can use the derivatives from the previous time steps
 * in our multistep method.
 */
void swap_buffers()
{
    double *tmp;

    tmp = dh2;
    dh2 = dh1;
    dh1 = dh;
    dh = tmp;

    tmp = du2;
    du2 = du1;
    du1 = du;
    du = tmp;

    tmp = dv2;
    dv2 = dv1;
    dv1 = dv;
    dv = tmp;
}

int t = 0;

/**
 * This is your step function!
 */
void step()
{
    // auto start = std::chrono::high_resolution_clock::now();
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> elapsed;

    // // First, we compute our ghost cells as we need them for our derivatives
    // start = std::chrono::high_resolution_clock::now();
    compute_ghost_horizontal();
    // end = std::chrono::high_resolution_clock::now();
    // elapsed = end - start;
    // std::cout << "compute_ghost_horizontal took " << elapsed.count() << " ms\n";

    // start = std::chrono::high_resolution_clock::now();
    compute_ghost_vertical();
    // end = std::chrono::high_resolution_clock::now();
    // elapsed = end - start;
    // std::cout << "compute_ghost_vertical took " << elapsed.count() << " ms\n";

    // // Next, we compute the derivatives of our fields
    // start = std::chrono::high_resolution_clock::now();
    compute_dh();
    // end = std::chrono::high_resolution_clock::now();
    // elapsed = end - start;
    // std::cout << "compute_dh took " << elapsed.count() << " ms\n";

    // start = std::chrono::high_resolution_clock::now();
    compute_du();
    // end = std::chrono::high_resolution_clock::now();
    // elapsed = end - start;
    // std::cout << "compute_du took " << elapsed.count() << " ms\n";

    // start = std::chrono::high_resolution_clock::now();
    compute_dv();
    // end = std::chrono::high_resolution_clock::now();
    // elapsed = end - start;
    // std::cout << "compute_dv took " << elapsed.count() << " ms\n";

    // We set the coefficients for our multistep method
    double a1 = 0.0, a2 = 0.0, a3 = 0.0;

    if (t == 0)
    {
        a1 = 1.0;
    }
    else if (t == 1)
    {
        a1 = 3.0 / 2.0;
        a2 = -1.0 / 2.0;
    }
    else
    {
        a1 = 23.0 / 12.0;
        a2 = -16.0 / 12.0;
        a3 = 5.0 / 12.0;
    }

    // // Finally, we compute the next time step using our multistep method
    // start = std::chrono::high_resolution_clock::now();
    multistep(a1, a2, a3);
    // end = std::chrono::high_resolution_clock::now();
    // elapsed = end - start;
    // std::cout << "multistep took " << elapsed.count() << " ms\n";

    // // We compute the boundaries for our fields
    // start = std::chrono::high_resolution_clock::now();
    compute_boundaries_horizontal();
    // end = std::chrono::high_resolution_clock::now();
    // elapsed = end - start;
    // std::cout << "compute_boundaries_horizontal took " << elapsed.count() << " ms\n";

    // start = std::chrono::high_resolution_clock::now();
    compute_boundaries_vertical();
    // end = std::chrono::high_resolution_clock::now();
    // elapsed = end - start;
    // std::cout << "compute_boundaries_vertical took " << elapsed.count() << " ms\n";

    // // We swap the buffers for our derivatives
    // start = std::chrono::high_resolution_clock::now();
    swap_buffers();
    // end = std::chrono::high_resolution_clock::now();
    // elapsed = end - start;
    // std::cout << "swap_buffers took " << elapsed.count() << " ms\n";

    t++;
}

/**
 * This is your transfer function! Since everything is running on the same node,
 * you don't need to do anything here.
 */
void transfer(double *h)
{
    return;
}

/**
 * This is your finalization function! Free whatever memory you've allocated.
 */
void free_memory()
{
    // TODO: Your code here

    free(dh);
    free(du);
    free(dv);

    free(dh1);
    free(du1);
    free(dv1);

    free(dh2);
    free(du2);
    free(dv2);
}
