#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include<mpi.h>
#include<pthread.h>

#include "grid.hpp"

// File I/O
// Liyu
// Output format is binary, to dump the binary file into ASCII format run: (replace N with n*ySize)
// hexdump -v -e 'N/4 "%10d "' -e '"\n"' outputfile > outputfile_ascii
void outputBGQ(GRID &localGrid, char* filename) { 
    MPI_File fh;
    // delete initial file 
    // http://stackoverflow.com/questions/25481661/how-to-replace-an-existing-file-in-mpi-with-mpi-file-open
    int err = MPI_File_open(MPI_COMM_WORLD, filename, 
    	MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if (err != MPI_SUCCESS)  {
        if (localGrid.myRank == 0)
            MPI_File_delete(filename, MPI_INFO_NULL);
        MPI_File_open(MPI_COMM_WORLD, filename, 
        	MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    }


    // MPI_File_open(MPI_COMM_WORLD, "outputfile", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    // int n = (int) sqrt((double) localGrid.totalRank);
    // MPI_Offset o = localGrid.xSize * n * localGrid.myOrigin[0] + localGrid.myOrigin[1];
    // Offset is: 
    // for (int i = 0; i < localGrid.ySize; i++){
    //     MPI_Offset offset = (localGrid.myOrigin[0] + i) * localGrid.xSize*((int)sqrt(localGrid.totalRank)) + localGrid.myOrigin[1];
    //     int* buf = localGrid.myUniverse[i]; // this
    //     int count = localGrid.xSize; // and this
    //     // try to write contigious memory to file
    //     MPI_Status* status;

    //     MPI_File_write_at(fh, offset*sizeof(MPI_INT), (void *) buf, count, MPI_INT, status);
    //     // printf("Rank %d (%d, %d), Offset %lld, Length %d\n", localGrid.myRank, localGrid.myOrigin[0], localGrid.myOrigin[1], offset, count);
    // }   
    // MPI_File_close(&fh);

	for (int i = 0; i < localGrid.xSize; i++){
	    MPI_Offset offset = (localGrid.myOrigin[0] + i) * localGrid.ySize*((int)sqrt(localGrid.totalRank)) + localGrid.myOrigin[1];
	    int* buf = localGrid.myUniverse[i]; // this
	    int count = localGrid.ySize; // and this
	    // try to write contigious memory to file
	    MPI_Status status;

	    MPI_File_write_at_all(fh, offset*sizeof(MPI_INT), (void *) buf, count, MPI_INT, &status);
	    // printf("Rank %d (%d, %d), Offset %lld, Length %d\n", localGrid.myRank, localGrid.myOrigin[0], localGrid.myOrigin[1], offset, count);
	}   
	MPI_File_close(&fh);

}