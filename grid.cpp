#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include<mpi.h>
#include<pthread.h>

#include "grid.hpp"

/*
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
    void swapGhost();
    void update();
    int mod (int a, int b);
};
*/


// ---------------------------------------------
// A method to initialize the grid and allocate moemories
void GRID::initialize(int size1, int size2, int mpiMyRank, int mpiTotalRank){
    
    myRank = mpiMyRank;
    totalRank = mpiTotalRank;
    int n = (int) sqrt((double) totalRank);
    xSize = size1 / n;
    ySize = size2 / n;

    
    // Initialize the origin of my universe (the global cordinate of myUniverse[0][0])
    myOrigin[0] = xSize * (myRank / n);
    myOrigin[1] = ySize * (myRank % n);
    
    // Find the neighbor ranks
    // 0 1   2   .... n-1
    // n n+1 n+2 .... 2n-1
    // ......
    // ......
    // n*(n-1) ....   n*n-1
    upperRank = (myRank < n) ? ((n-1)*n + myRank) : (myRank - n);
    lowerRank = (myRank >= (n-1)*n) ? (myRank - (n-1)*n) : (myRank + n);
    leftRank = (myRank % n == 0) ? (myRank + n - 1) : (myRank - 1);
    rightRank = (myRank % n == n - 1) ? (myRank - n + 1) : (myRank + 1);
    
    upperLeftRank = (upperRank % n == 0) ? (upperRank + n - 1) : (upperRank - 1);
    upperRightRank = (upperRank % n == n - 1) ? (upperRank - n + 1) : (upperRank + 1);
    lowerLeftRank = (lowerRank % n == 0) ? (lowerRank + n - 1) : (lowerRank - 1);
    lowerRightRank = (lowerRank % n == n - 1) ? (lowerRank - n + 1) : (lowerRank + 1);

    //if (myRank == 0)
        //printf("Rank %d : upper %d, lower %d, left %d, right %d\n", myRank, upperRank, lowerRank, leftRank, rightRank);
    
    // Allocate contiguous memory for myUniverse
    arrayHead = (int *) malloc(xSize*ySize*sizeof(int)); 
    myUniverse = (int**) malloc(xSize*sizeof(int*));
    
    // Assign head-of-row pointers to myUniverse
    for (int i = 0; i < xSize; i++){
        myUniverse[i] = &arrayHead[ySize*i];
    }
    // Allocate memory for ghost arrays
    upperGhost = (int *) malloc(ySize * sizeof(int));
    lowerGhost = (int *) malloc(ySize * sizeof(int));
    leftGhost = (int *) malloc(xSize * sizeof(int));
    rightGhost = (int *) malloc(xSize * sizeof(int));
}

void GRID::freeMemory(){
    free(arrayHead);
    free(upperGhost);
    free(lowerGhost);
    free(leftGhost);
    free(rightGhost);
}

// Converting local coordinate to global coordinate    
void GRID::local2Global(int* localCord, int* globalCord){
    globalCord[0] = localCord[0] + myOrigin[0];
    globalCord[1] = localCord[1] + myOrigin[1];
}

// TODO: CJ
void GRID::swapGhost(int tag){
    MPI_Datatype col_type;
    MPI_Type_vector(xSize, 1, ySize, MPI_INT, &col_type);
    MPI_Type_commit(&col_type);
    
    MPI_Request sendreq;
    MPI_Status status;
    MPI_Isend(myUniverse[0], ySize, MPI_INT, upperRank, tag, MPI_COMM_WORLD, &sendreq);
    MPI_Isend(myUniverse[xSize-1], ySize, MPI_INT, lowerRank, tag, MPI_COMM_WORLD, &sendreq);
    MPI_Isend(&myUniverse[0][0], 1, col_type, leftRank, tag, MPI_COMM_WORLD, &sendreq);
    MPI_Isend(&myUniverse[0][ySize-1], 1, col_type, rightRank, tag, MPI_COMM_WORLD, &sendreq);
    
    // four corners
    MPI_Isend(&myUniverse[0][0], 1, MPI_INT, upperLeftRank, tag, MPI_COMM_WORLD, &sendreq);
    MPI_Isend(&myUniverse[0][ySize-1], 1, MPI_INT, upperRightRank, tag, MPI_COMM_WORLD, &sendreq);
    MPI_Isend(&myUniverse[xSize-1][0], 1, MPI_INT, lowerLeftRank, tag, MPI_COMM_WORLD, &sendreq);
    MPI_Isend(&myUniverse[xSize-1][ySize-1], 1, MPI_INT, lowerRightRank, tag, MPI_COMM_WORLD, &sendreq);



    MPI_Request recvreq;
    MPI_Irecv(lowerGhost, ySize, MPI_INT, lowerRank, tag, MPI_COMM_WORLD, &recvreq);
    MPI_Wait(&recvreq, &status);
    MPI_Irecv(upperGhost, ySize, MPI_INT, upperRank, tag, MPI_COMM_WORLD, &recvreq);
    MPI_Wait(&recvreq, &status);
    MPI_Irecv(leftGhost, xSize, MPI_INT, leftRank, tag, MPI_COMM_WORLD, &recvreq);
    MPI_Wait(&recvreq, &status);
    MPI_Irecv(rightGhost, xSize, MPI_INT, rightRank, tag, MPI_COMM_WORLD, &recvreq);
    MPI_Wait(&recvreq, &status);

    MPI_Irecv(&upperLeftCorner, 1, MPI_INT, upperLeftRank, tag, MPI_COMM_WORLD, &recvreq);
    MPI_Wait(&recvreq, &status);
    MPI_Irecv(&upperRightCorner, 1, MPI_INT, upperRightRank, tag, MPI_COMM_WORLD, &recvreq);
    MPI_Wait(&recvreq, &status);
    MPI_Irecv(&lowerLeftCorner, 1, MPI_INT, lowerLeftRank, tag, MPI_COMM_WORLD, &recvreq);
    MPI_Wait(&recvreq, &status);
    MPI_Irecv(&lowerRightCorner, 1, MPI_INT, lowerRightRank, tag, MPI_COMM_WORLD, &recvreq);
    MPI_Wait(&recvreq, &status);
}


// TODO: Frank
void GRID::update() // only a single tick is enough
{
    /*
    Pseudocode:
    For all the lattice points{
        If (at Grain Boundary) {
            - Randomly select a neighbor grain ID
            - Compute energy difference: Original Energy - Replaced Energy ==> \Delta E
            - Make decision: accept new grain ID or not?
                    \Delta E < 0 => accept
            - Modify grain ID
        
        }
    }
    */
    //printf("xSize: %d",xSize);
    //printf("ySize: %d",ySize);
    //loop through all the rows
    for (int i =0; i < xSize; i++)
    {
        // iterate all columns for current row
        for (int j=0; j<ySize;j++)
        {
            int map[8][2];
            for (int atK = 0; atK < 8; atK++)
            {
                map[atK][0] = 0;
                map[atK][1] = 0;
                //printf("atK: %d with 0 value: %d 1 value: %d",atK,map[atK][0],map[atK][1]);
            }
            int numberOfGrains = 0;
            int currentEnergy = 0;
             /* count if current cell is on grain boundary */
            int i_u = mod(i-1, xSize);//upper row
            //printf("i_u is: %d",i_u);
            int i_d = mod(i+1, xSize);//lower row
            //printf("i_d is: %d",i_d);
            int j_r = mod(j+1, ySize);//right column
            //printf("j_r is: %d",j_r);
            int j_l = mod(j-1, ySize);//left column
            //printf("j_l is: %d",j_l);
            // 0 - upper left 1 - upper 2 - upper right
            // 3 - left                 4- right
            // 5 - lower left 6 - lower 7 - lower right
            int neighbourGrains[8];
            neighbourGrains[0] = myUniverse[i_u][j_l];
            neighbourGrains[1] = myUniverse[i_u][j];
            neighbourGrains[2] = myUniverse[i_u][j_r];
            neighbourGrains[3] = myUniverse[i][j_l];
            neighbourGrains[4] = myUniverse[i][j_r];
            neighbourGrains[5] = myUniverse[i_d][j_l];
            neighbourGrains[6] = myUniverse[i_d][j];
            neighbourGrains[7] = myUniverse[i_d][j_r];

            //specical consideration for upper row
            // the first row
            if (i == 0)                
            {
                neighbourGrains[0] = upperGhost[j_l];
                neighbourGrains[1] = upperGhost[j];
                neighbourGrains[2] = upperGhost[j_r];
            }

            //the last row
            if (i == xSize -1)                
            {
                neighbourGrains[5] = lowerGhost[j_l];
                neighbourGrains[6] = lowerGhost[j];
                neighbourGrains[7] = lowerGhost[j_r];
            }

            //the first column
            if (j == 0)
            {
                neighbourGrains[0] = leftGhost[i_u];
                neighbourGrains[3] = leftGhost[i];
                neighbourGrains[5] = leftGhost[i_d];
            }
            
            //the last column
            if (j == ySize -1)                
            {
                neighbourGrains[2] = rightGhost[i_u];
                neighbourGrains[4] = rightGhost[i];
                neighbourGrains[7] = rightGhost[i_d];
            }
            
            //if upper left corner
            if ((i == 0) && (j==0))
            neighbourGrains[0] = upperLeftCorner;
            
            //if upper right corner
            if ((i == 0) && (j== ySize-1))
            neighbourGrains[2] = upperRightCorner;
            
            //if lower left corner
            if ((i==xSize-1) && j==0)
            neighbourGrains[5] = lowerLeftCorner;
            
            //if lower right corner
            if (i==xSize-1 && j== ySize-1)
            neighbourGrains[7] = lowerRightCorner;

            for (int kk = 0; kk<8; kk++)
            {
                if (neighbourGrains[kk] != myUniverse[i][j])
                {
                    bool found = false;
                    for (int atK = 0; ((atK < 8) && (found == false)); atK++)
                    {
                        if (map[atK][0] == neighbourGrains[kk])
                        {
                            map[atK][1]++;
                            found = true;
                        }
                        if (map[atK][0] == 0)
                        {
                            map[atK][0] = neighbourGrains[kk];
                            map[atK][1]++;
                            numberOfGrains ++;
                            found = true;
                        }
                    }
                }  
                else
                    currentEnergy++;
            }

            // grain boundary true
            if (numberOfGrains > 0)
            {
                //local2Global(&i,&j);
                //double r = GenVal(i + myOrigin[0] + j + myOrigin[1]);
                double r = (double)rand()/(double)RAND_MAX;
                // randomly select a neighbour grain
                double range = 1.0/(numberOfGrains);
                bool replaced = false;
                for (int atK = 0; ((atK < numberOfGrains)&&(replaced==false)); atK++)
                {
                    if (r <= (atK +1) * range)
                    {
                        //compute delta E
                        int deltaE = map[atK][1] - currentEnergy;
                        //if energy difference < 0, then replace current grain ID with neighbour grain ID
                        if (deltaE < 0)
                        {
                        myUniverse[i][j] = map[atK][0];
                        }
                        replaced = true;
                    }
                }
            }
        }
    }
}


void GRID::updateCJ() // only a single tick is enough
{
    /*
    Pseudocode:
    For all the lattice points{
        If (at Grain Boundary) {
            - Randomly select a neighbor grain ID
            - Compute energy difference: Original Energy - Replaced Energy ==> \Delta E
            - Make decision: accept new grain ID or not?
                    \Delta E < 0 => accept
            - Modify grain ID
        
        }
    }
    */
    int iter = 0;
    srand(myRank + (int) time(0));

    while(iter < xSize * ySize){
        int x = rand() % xSize;
        int y = rand() % ySize;
        // int x = floor(GenVal(myRank) * xSize);
        // int y = floor(GenVal(myRank) * ySize);
        int upper, lower, left, right, upperLeft, upperRight, lowerLeft, lowerRight;

        // find the 8 neighbors !!
        upper = (x == 0) ? upperGhost[y]: myUniverse[x - 1][y];
        lower = (x == xSize-1) ? lowerGhost[y]: myUniverse[x + 1][y];
        left = (y == 0) ? leftGhost[x]: myUniverse[x][y - 1];
        right = (y == ySize-1) ? rightGhost[x]: myUniverse[x][y + 1];

        if (x == 0 && y == 0)
            upperLeft = upperLeftCorner;
        else if (x == 0) // first row
            upperLeft = upperGhost[y-1];
        else if  (y == 0) // first col
            upperLeft = leftGhost[x-1];
        else 
            upperLeft = myUniverse[x - 1][y - 1];

        if (x == 0 && y == ySize - 1)
            upperRight = upperRightCorner;
        else if (x == 0) // first row
            upperRight = upperGhost[y + 1];
        else if  (y == ySize - 1) // last col
            upperRight = rightGhost[x - 1];
        else 
            upperRight = myUniverse[x - 1][y + 1];

        if (x == xSize-1 && y == 0)
            lowerLeft = lowerLeftCorner;
        else if (x == xSize-1) // last row
            lowerLeft = lowerGhost[y - 1];
        else if  (y == 0) // first col
            lowerLeft = leftGhost[x + 1];
        else 
            lowerLeft = myUniverse[x + 1][y - 1];

        if (x == xSize-1 && y == ySize - 1)
            lowerRight = lowerRightCorner;
        else if (x == xSize-1) // last row
            lowerRight = lowerGhost[y + 1];
        else if  (y == ySize-1) // last col
            lowerRight = rightGhost[x + 1];
        else 
            lowerRight = myUniverse[x + 1][y + 1];

        // count how many grains are different
        int neib[8] = {upper, lower, left, right, upperLeft, upperRight, lowerLeft, lowerRight};
        int myID = myUniverse[x][y];
        bool GB = false;
        int energy = 0;
        int diffCount = 0;
        for (int i = 0; i < 8; i++){
            if (neib[i] != myID){
                GB = true;
                energy ++;
            }
        }
        if (GB == true){
            //printf("I'm at GB!\n");
            int newID = neib[rand() % 8];
            while (newID == myID){
                newID = neib[rand() % 8];
            }
            int newEnergy = 0;
            for (int i = 0; i < 8; i++){
                if (neib[i] != newID){
                    newEnergy ++;
                }
            }
            int kT = 1.3806488e-23*423;
            double r = double(rand())/double(RAND_MAX);
            if (energy - newEnergy >=0)
                myUniverse[x][y] = newID;
            else if ( r < exp( ((double) (energy - newEnergy))/kT ) ) {
                myUniverse[x][y] = newID;
                //printf("I fliped from %d to %d on (%d, %d) !\n", myID, newID, x, y);
            }
        }
        iter ++;
    }
}

pthread_mutex_t **locks;
pthread_mutex_t inputlock;
pthread_barrier_t barrier;

void GRID::threadUpdate(int nT){
    pthread_barrier_init(&barrier, NULL, nT);
    locks = (pthread_mutex_t **) malloc((xSize+2)*sizeof( pthread_mutex_t * ));
    for (int m = 0; m < xSize + 2; m++){
        locks[m] = (pthread_mutex_t *) malloc((ySize+2)*sizeof( pthread_mutex_t ));
        for (int n = 0; n < ySize + 2; n++){
            pthread_mutex_init( &locks[m][n], NULL );
        }
    }
    pthread_mutex_init(&inputlock, NULL);
    pthread_t tid[nT];
    //printf("Rank %d: what the fuck am I doing?", myRank);

    for (int i = 0; i < nT; i++){
        TInput* input = (TInput *) malloc(sizeof(TInput));
        input->tid = i;
        input->totalThread = nT;
        input->myUniverse = myUniverse;
        input->upperGhost = upperGhost;
        input->lowerGhost = lowerGhost;
        input->leftGhost = leftGhost;
        input->rightGhost = rightGhost;

        input->upperLeftCorner = upperLeftCorner;
        input->upperRightCorner = upperRightCorner;
        input->lowerLeftCorner = lowerLeftCorner;
        input->lowerRightCorner = lowerRightCorner;

        input->xSize = xSize;
        input->ySize = ySize;
        input->myRank = myRank;
        input->totalRank = totalRank;
        //printf("Rank %d: trying to create thread %d\n", myRank, i);
        int rc = pthread_create(&tid[i], NULL, threadHelper, (void *) input);
        //printf("Rank %d: Successfully created thread %d \n", myRank, i);
    }
    for (int j = 0; j < nT; j++){
        pthread_join(tid[j], NULL);
            //if (mpi_myrank==0) printf("Successfully joined TID: %d \n",tid[i]);
    }
}



void* threadHelper(void *arg){
    TInput* input = (TInput *) arg;

    pthread_mutex_lock(&inputlock);
    int tid = input->tid;
    int totalThread = input->totalThread;
    int** myUniverse = input->myUniverse;
    int* upperGhost = input->upperGhost;
    int* lowerGhost = input->lowerGhost;
    int* leftGhost = input->leftGhost;
    int* rightGhost = input->rightGhost;

    int upperLeftCorner = input->upperLeftCorner;
    int upperRightCorner = input->upperRightCorner;
    int lowerLeftCorner = input->lowerLeftCorner;
    int lowerRightCorner = input->lowerRightCorner;

    int xSize = input->xSize;
    int ySize = input->ySize;
    int myRank = input->myRank;
    int totalRank = input->totalRank;
    free(input);
    pthread_mutex_unlock(&inputlock);

    //if (myRank ==0) printf("TID %d: I'm in!\n", tid);
    int iter = 0;
    srand(myRank*tid + tid + (int) time(0));
    
    int rowsPerThread = xSize / totalThread;

    while(iter < rowsPerThread * ySize){
        int x = rand() % rowsPerThread +  tid*rowsPerThread;
        int y = rand() % ySize;
        //if (myRank ==0) 
            //printf("Rank%d: TID %d: The %dth iteration, (x, y) = (%d, %d)\n", myRank, tid, iter, x, y);
        // int x = floor(GenVal(myRank) * xSize);
        // int y = floor(GenVal(myRank) * ySize);
        int upper, lower, left, right, upperLeft, upperRight, lowerLeft, lowerRight;

        for (int mm = x; mm < x+3; mm++)
            for (int nn = y; nn < y+3; nn++)
                pthread_mutex_lock(&locks[mm][nn]);

        // find the 8 neighbors !!
        upper = (x == 0) ? upperGhost[y]: myUniverse[x - 1][y];
        lower = (x == xSize-1) ? lowerGhost[y]: myUniverse[x + 1][y];
        left = (y == 0) ? leftGhost[x]: myUniverse[x][y - 1];
        right = (y == ySize-1) ? rightGhost[x]: myUniverse[x][y + 1];

        if (x == 0 && y == 0)
            upperLeft = upperLeftCorner;
        else if (x == 0) // first row
            upperLeft = upperGhost[y-1];
        else if  (y == 0) // first col
            upperLeft = leftGhost[x-1];
        else 
            upperLeft = myUniverse[x - 1][y - 1];

        if (x == 0 && y == ySize - 1)
            upperRight = upperRightCorner;
        else if (x == 0) // first row
            upperRight = upperGhost[y + 1];
        else if  (y == ySize - 1) // last col
            upperRight = rightGhost[x - 1];
        else 
            upperRight = myUniverse[x - 1][y + 1];

        if (x == xSize-1 && y == 0)
            lowerLeft = lowerLeftCorner;
        else if (x == xSize-1) // last row
            lowerLeft = lowerGhost[y - 1];
        else if  (y == 0) // first col
            lowerLeft = leftGhost[x + 1];
        else 
            lowerLeft = myUniverse[x + 1][y - 1];

        if (x == xSize-1 && y == ySize - 1)
            lowerRight = lowerRightCorner;
        else if (x == xSize-1) // last row
            lowerRight = lowerGhost[y + 1];
        else if  (y == ySize-1) // last col
            lowerRight = rightGhost[x + 1];
        else 
            lowerRight = myUniverse[x + 1][y + 1];

        // count how many grains are different
        int myID = myUniverse[x][y];
        for (int mm = x; mm < x+3; mm++)
            for (int nn = y; nn < y+3; nn++)
                pthread_mutex_unlock(&locks[mm][nn]);
        
        int neib[8] = {upper, lower, left, right, upperLeft, upperRight, lowerLeft, lowerRight};
        bool GB = false;
        int energy = 0;
        int diffCount = 0;
        for (int i = 0; i < 8; i++){
            if (neib[i] != myID){
                GB = true;
                energy ++;
            }
        }
        if (GB == true){
            //printf("I'm at GB!\n");
            int newID = neib[rand() % 8];
            while (newID == myID){
                newID = neib[rand() % 8];
            }
            int newEnergy = 0;
            for (int i = 0; i < 8; i++){
                if (neib[i] != newID){
                    newEnergy ++;
                }
            }
            int kT = 1.3806488e-23*423;
            double r = double(rand())/double(RAND_MAX);
            if (energy - newEnergy >=0){
                pthread_mutex_lock(&locks[x+1][y+1]);
                myUniverse[x][y] = newID;
                pthread_mutex_unlock(&locks[x+1][y+1]);
            }
            else if ( r < exp( ((double) (energy - newEnergy))/kT ) ) {
                pthread_mutex_lock(&locks[x+1][y+1]);
                myUniverse[x][y] = newID;
                pthread_mutex_unlock(&locks[x+1][y+1]);
                //printf("I fliped from %d to %d on (%d, %d) !\n", myID, newID, x, y);
            }
        }
        iter ++;
        pthread_barrier_wait(&barrier);
        //if (myRank ==0) printf("TID %d: starting the %dth iteration\n", tid, iter);
    }
    pthread_exit(0);
}

/* mod function */
int mod (int a, int b)
{
    if (a < 0)
        return (b-1);
    else
        if (a >= b)
            return 0;
        else
            return a;
}


