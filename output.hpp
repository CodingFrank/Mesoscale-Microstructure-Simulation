# ifndef OUTPUT_HPP
# define OUTPUT_HPP

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include<mpi.h>
#include<pthread.h>

#include "grid.hpp"

void outputBGQ(GRID &localGrid, char* );

# endif