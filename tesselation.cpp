#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include<mpi.h>
#include<pthread.h>
#include <cstdlib>
#include <ctime>

#include "grid.hpp"

// Tesselation

typedef struct {
    int x; 
    int y;
} cord;

/*
    Initialize the domain by Voronoi Tessellation
    Written by CJ

    Pseudocode:

    tesselation(GRID &grid, int numSeeds){
        if (rank==0)
            Randomly generate N seeds in the global domain \
                store them in the array seeds[N][2];
        Broadcast all N seeds from rank 0 to all the ranks;
        
        for all the lattice points in my local grid{
            int minDistance = INT_MAX;
            int minID = 0;
            for (int i=0; i<N; i++){
                Compute the Distance from seed[i] to the current lattice point;
                if (Distance < minDistance){
                    minDistance = Distance;
                    minID = i;
                }
            }
            Modify the grain ID in current lattice point as minID;
        }
    }    

*/



void tesselation(GRID &grid, int numSeeds){

    while (numSeeds % grid.totalRank != 0) ++numSeeds;

    int numMySeeds = numSeeds / grid.totalRank;

    // cord* mySeeds = new cord[numMySeeds];
    // cord* allSeeds = new cord[numSeeds];

    int cordx[numSeeds];
    int cordy[numSeeds];

    srand((unsigned)time(0));

    for (int i = 0; i < numSeeds; i++){
        cordx[i] = rand() % (grid.xSize * (int) sqrt(grid.totalRank));
        cordy[i] = rand() % (grid.ySize * (int) sqrt(grid.totalRank));
    }

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast( cordx, numSeeds, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast( cordy, numSeeds, MPI_INT, 0, MPI_COMM_WORLD );

    //MPI_Allgather(mySeeds, numMySeeds, MPI_INT, allSeeds, numSeeds, MPI_INT, MPI_COMM_WORLD);
    
    if (grid.myRank == 0) printf("Tesselation is starting!\n");

    for (int x = 0; x < grid.xSize; x++){
        for (int y = 0; y < grid.ySize; y++){
            int local[2] = {x,y};
            int global[2];
            grid.local2Global(local, global);
            int minDistance = pow(global[0] - cordx[0], 2) 
                        + pow(global[1] - cordy[0] , 2);
            int minID = 1;
            for (int i = 0; i < numSeeds; i++){
                int dist = pow(global[0] - cordx[i], 2) 
                        + pow(global[1] - cordy[i] , 2);
                if (dist < minDistance){
                    minDistance = dist;
                    minID = i + 1;
                }
            }
            grid.myUniverse[x][y] = minID;
        }
    }
    
    //delete[] mySeeds;
    //delete[] allSeeds;
}