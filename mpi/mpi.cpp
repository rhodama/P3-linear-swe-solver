#include <mpi.h>
#include <stdlib.h>
#include <cstdio>
#include <cstring> 

#include "../common/common.hpp"
#include "../common/solver.hpp"
int nx_all, ny_all;
int nx,ny;
int rank,num_procs;

double *h, *u, *v, *dh, *du, *dv, *dh1, *du1, *dv1, *dh2, *du2, *dv2;

double H, g, dx, dy, dt;
int t=0;
/**
 * This is your initialization function! It is very similar to the one in
 * serial.cpp, but with some difference. Firstly, only the process with rank 0
 * is going to actually generate the initial conditions h0, u0, and v0, so all
 * other processes are going to get nullptrs. Therefore, you'll need to find some
 * way to scatter the initial conditions to all processes. Secondly, now the
 * rank and num_procs arguments are passed to the function, so you can use them
 * to determine which rank the node running this process has, and how many
 * processes are running in total. This is useful to determine which part of the
 * domain each process is going to be responsible for.
 */
void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    // TODO: Your code here

    rank = rank_;
    num_procs =num_procs_;
    int size;
   

    nx_all = nx_;
    ny_all = ny_;
    nx = nx_all / num_procs;
    if (rank == num_procs - 1) {
    nx += nx_all % num_procs;  // The last process handles the remaining rows
    }
    ny =ny_all;
    size=nx *ny;
    // We allocate memory for the derivatives
    h = (double *)calloc((nx + 1) * (ny + 1), sizeof(double));
    if (h == NULL) {
        fprintf(stderr, "Failed to allocate memory for h\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    u = (double *)calloc((nx + 1) * ny, sizeof(double));
    if (u == NULL) {
        fprintf(stderr, "Failed to allocate memory for u\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    v = (double *)calloc(nx * (ny + 1), sizeof(double));
    if (v == NULL) {
        fprintf(stderr, "Failed to allocate memory for v\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    dh = (double *)calloc(size, sizeof(double));
    du = (double *)calloc(size, sizeof(double));
    dv = (double *)calloc(size, sizeof(double));

    dh1 = (double *)calloc(size, sizeof(double));
    du1 = (double *)calloc(size, sizeof(double));
    dv1 = (double *)calloc(size, sizeof(double));

    dh2 = (double *)calloc(size, sizeof(double));
    du2 = (double *)calloc(size, sizeof(double));
    dv2 = (double *)calloc(size, sizeof(double));

    H = H_;
    g = g_;

    dx = length_ / nx_all;
    dy = width_ / ny_all;

    dt = dt_;

    
    int *sendcounts_h = (int*)malloc(num_procs * sizeof(int));
    int *displs_h = (int*)malloc(num_procs * sizeof(int));

    int *sendcounts_u = (int*)malloc(num_procs * sizeof(int));
    int *displs_u = (int*)malloc(num_procs * sizeof(int));

    int *sendcounts_v = (int*)malloc(num_procs * sizeof(int));
    int *displs_v = (int*)malloc(num_procs * sizeof(int));    

    displs_u[0] = 0;
    displs_h[0] =0;
    displs_v[0] =0;

    for(int i = 0; i < num_procs; i++) {
        sendcounts_h[i] = (ny+1) * (nx_all / num_procs);

        sendcounts_u[i] = (ny) * ((nx_all) / num_procs);

        sendcounts_v[i] = (ny+1) * ((nx_all) / num_procs);

        if(i == num_procs - 1) {
            sendcounts_h[i] += (ny+1) * (nx_all % num_procs+1);
            sendcounts_u[i] += (ny) * ((nx_all) % num_procs+1);
            sendcounts_v[i] += (ny+1) * ((nx_all) % num_procs);
        }
        if(i > 0) {
            displs_h[i] = displs_h[i-1] + sendcounts_h[i-1];
            displs_u[i] = displs_u[i-1] + sendcounts_u[i-1];
            displs_v[i] = displs_v[i-1] + sendcounts_v[i-1];
        }
    }
    
    
    MPI_Scatterv(h0, sendcounts_h, displs_h, MPI_DOUBLE,
                h, (nx + 1) * (ny + 1), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
    MPI_Scatterv(u0, sendcounts_u, displs_u, MPI_DOUBLE,
                u, (nx + 1) * ny, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
    MPI_Scatterv(v0, sendcounts_v, displs_v, MPI_DOUBLE,
                v, nx * (ny + 1), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
 
    free(sendcounts_h);
    free(sendcounts_u);
    free(sendcounts_v);

    free(displs_h);
    free(displs_u);
    free(displs_v);

}

void exchange_u_h_cells()
{
    MPI_Status status;
    int prev_rank = (rank - 1 + num_procs) % num_procs;
    int next_rank = (rank + 1) % num_procs;
    
   if((rank!=0)&&(rank!=num_procs-1)){

        MPI_Sendrecv(
            &u(0, 0), ny, MPI_DOUBLE, prev_rank, 0,
            &u(nx,0), ny, MPI_DOUBLE, next_rank, 0,
            MPI_COMM_WORLD, &status
        );  
         MPI_Sendrecv(
            &h(0, 0), ny, MPI_DOUBLE, prev_rank, 0,
            &h(nx,0), ny, MPI_DOUBLE, next_rank, 0,
            MPI_COMM_WORLD, &status
        );  
   }
   else if(rank==0){

        MPI_Recv(
            &u(nx,0), ny, MPI_DOUBLE,1,0, MPI_COMM_WORLD, &status
        );    
         MPI_Recv(
            &h(nx,0), ny, MPI_DOUBLE,1,0, MPI_COMM_WORLD, &status
        ); 
   }
   else if(rank ==num_procs-1)
   {

        MPI_Send(
            &u(0,0), ny, MPI_DOUBLE, num_procs-2, 0, MPI_COMM_WORLD
        );
        MPI_Send(
            &h(0,0), ny, MPI_DOUBLE, num_procs-2, 0, MPI_COMM_WORLD
        );

   }

}

void exchange_ghost_cells()
{
    

    for (int i = 0; i < nx; i++)
    {
        h(i, ny) = h(i, 0);
    }
     for (int j = 0; j < ny; j++)
    {
        h(nx, j) = h(0, j);
    }
    for (int j = 0; j < ny; j++)
    {
        u(0, j) = u(nx, j);
    }
    for (int i = 0; i < nx; i++)
    {
        v(i, 0) = v(i, ny);
    }

}

void compute_dh_du_dv()
{
    
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
            du(i, j) = -g * dh_dx(i, j);
            dv(i, j) = -g * dh_dy(i, j);
        }
    }
}

void multistep(double a1, double a2, double a3)
{
    
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
            u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
            v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
        }
    }
}

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
/**
 * This is your step function! It is very similar to the one in serial.cpp, but
 * now the domain is divided among the processes, so you'll need to find some
 * way to communicate the ghost cells between processes.
 */
void step()
{
    // TODO: Your code here
    // First, we compute our ghost cells as we need them for our derivatives
    exchange_ghost_cells();
    exchange_u_h_cells();

    // Next, we compute the derivatives of our fields
    compute_dh_du_dv();


    // We set the coefficients for our multistep method
    double a1, a2, a3;

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

    // Finally, we compute the next time step using our multistep method
    multistep(a1, a2, a3);


    swap_buffers();

    t++;
}

/**
 * This is your transfer function! Similar to what you did in gpu.cu, you'll
 * need to get the data from the computers you're working on (there it was
 * the GPU, now its a bunch of CPU nodes), and send them all back to the process
 * which is actually running the main function (then it was the CPU, not it's
 * the node with rank 0).
 */
void transfer(double *h_recv)
{
    // TOOD: Your code here
    int *recvcounts = (int*)malloc(num_procs * sizeof(int));
    int sendcounts_1; 

    int *displs = (int*)malloc(num_procs * sizeof(int));


nx = nx_all / num_procs;
if (rank == num_procs - 1) {
    nx += nx_all % num_procs;  // The last process handles the remaining rows
}
ny = ny_all;

displs[0] = 0;
// Calculate the data amount and displacement for each process
for (int i = 0; i < num_procs; i++) {
    // Calculate the number of rows each process is responsible for
    recvcounts[i] = (ny + 1) * (nx_all / num_procs);
    if (i == num_procs - 1) {
        recvcounts[i] += (ny + 1) * (nx_all % num_procs + 1);  // Add remaining rows to the last process
    }
    if (i > 0) {
        displs[i] = displs[i - 1] + recvcounts[i - 1];  // Calculate displacement based on previous process
    }
}


    if(rank!=num_procs-1){
        sendcounts_1 = nx*(ny+1);
    }
    else if(rank==num_procs-1){
        sendcounts_1 = (nx+1)*(ny+1);
    }

// Use MPI_Gatherv to collect data to rank 0
MPI_Gatherv(
    h,               
    sendcounts_1,   
    MPI_DOUBLE,      
    h_recv,          
    recvcounts,     
    displs,          
    MPI_DOUBLE,      
    0,               
    MPI_COMM_WORLD   
);

    
    delete[] recvcounts;
    delete[] displs;
    
}

/**
 * This is your finalization function! Since different nodes are going to be
 * initializing different chunks of memory, make sure to check which node
 * is running the code before you free some memory you haven't allocated, or
 * that you've actually freed memory that you have.
 */
void free_memory()
{
    // TODO: Your code here
// Release the derivative buffers allocated by all processes
        free(dh);

        free(du);
        free(dv);
    
        free(dh1);
        free(du1);
        free(dv1);
    
        free(dh2);
        free(du2);
        free(dv2);
        free(h);
        free(u);
        free(v);
}

