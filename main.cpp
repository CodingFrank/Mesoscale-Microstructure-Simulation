/***************************************************************************/
/* PPC Asssignment 56          **********************************************/
/* Chengjian Zheng            **********************************************/
/* RIN: 661-23-8526           **********************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include<mpi.h>
#include<pthread.h>

#include "grid.hpp"
#include "output.hpp"
#include "tesselation.hpp"
/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
    if (argc != 6)
    {
        printf("Invalid argument list. Use mc.out <Number of threads per rank> <Domain size> <Number of updates> <Output Interval> <Initial grian size>\n");
        return EXIT_FAILURE;
    }
    int totalThread = atoi(argv[1]);
    int Dsize = atoi(argv[2]);
    int Nstep = atoi(argv[3]);
    int writeInterval = atoi(argv[4]);
    double initialSize = atof(argv[5]);
    int initialNumGrains = (int) pow(((double) Dsize) / initialSize , 2);

    int mpi_myrank;
    int mpi_commsize;
    // Example MPI startup and using CLCG4 RNG
    MPI_Init(NULL, NULL);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
    double startTime;
    startTime = MPI_Wtime();

    GRID grid;
    grid.initialize(Dsize, Dsize, mpi_myrank, mpi_commsize);
    if (mpi_myrank == 0) 
        printf("Initialization done!\n");
    tesselation(grid, initialNumGrains);
    MPI_Barrier(MPI_COMM_WORLD);

    double timeNow;
    timeNow = MPI_Wtime();


    if (mpi_myrank == 0) 
        printf("Tessellation done! Time spent is %f s\n", timeNow - startTime);
    outputBGQ(grid, "initial.mcdat");
    //grid.update();



    for (int i = 0; i < Nstep; i++){
        if (mpi_myrank == 0) 
            printf("------------- The %d th step -----------\n", i);

        grid.swapGhost(0);
        if (mpi_myrank == 0) 
            printf("Swapping ghost done!\n");
        
        if (totalThread == 0)
            grid.updateCJ();
        else
            grid.threadUpdate(totalThread);

        MPI_Barrier(MPI_COMM_WORLD);
        if (mpi_myrank == 0)
            printf("Update done!\n");

        if (i % writeInterval == 0){
            char filename[20];
            sprintf(filename, "snapshot%d.mcdat",i);
            outputBGQ(grid, filename);
            if (mpi_myrank == 0)
                printf("Written output to file!\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (mpi_myrank == 0) 
        printf("Update done! Time spent is %f s\n", timeNow - startTime);
    
    outputBGQ(grid, "final.mcdat");

    grid.freeMemory();
    if (mpi_myrank == 0) printf("Freed the memory, call it a day!\n");
    MPI_Finalize();
    return 0;
}
