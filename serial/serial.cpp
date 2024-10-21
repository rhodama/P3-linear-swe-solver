#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../common/common.hpp"
#include "../common/solver.hpp"

/**
 * This is your initialization function!
 */
void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    // TODO: Your code here
}

/**
 * This is your step function!
 */
void step()
{
    // TODO: Your code here
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
}