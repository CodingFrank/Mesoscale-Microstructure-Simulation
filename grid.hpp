# ifndef GRID_HPP
# define GRID_HPP

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include<mpi.h>
#include<pthread.h>

class GRID {
    // ---------------------------------------------
    // variable definitions needs to be communicated
public:
    int xSize;
    int ySize;
    int myRank;
    int totalRank;
    
    int myOrigin[2];  
    int leftRank;
    int rightRank;
    int upperRank;
    int lowerRank;
    
    int upperLeftRank;
    int upperRightRank;
    int lowerLeftRank;
    int lowerRightRank;

    // all the universe is stored in arrayHead
    int* arrayHead;
    // myUniverse is an array of head-of-row pointers
    // Access the lattice point [i,j] by myUniverse[i][j]
    int** myUniverse;
    // the ghost rows and columns are stored in these four arrays
    int* upperGhost;
    int* lowerGhost;
    int* leftGhost;
    int* rightGhost;
    
    // the ghost corners
    int upperLeftCorner;
    int upperRightCorner;
    int lowerLeftCorner;
    int lowerRightCorner;
    
    void initialize(int size1, int size2, int mpiMyRank, int mpiTotalRank);
    void freeMemory();
    void local2Global(int* localCord, int* globalCord);
    void swapGhost(int);
    void update();
    void updateCJ();
    void threadUpdate(int nT);
};

typedef struct ThreadInput{
    int tid;
    int totalThread;
    
    int** myUniverse;
    int* upperGhost;
    int* lowerGhost;
    int* leftGhost;
    int* rightGhost;

    int upperLeftCorner;
    int upperRightCorner;
    int lowerLeftCorner;
    int lowerRightCorner;

    int xSize;
    int ySize;
    int myRank;
    int totalRank;
} TInput;

void* threadHelper(void *input);

int mod (int a, int b);
//void tesselation(GRID &localGrid, int seeds);

# endif
